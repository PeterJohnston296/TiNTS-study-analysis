#!/usr/bin/env python3
import argparse, csv, gzip, json, math, os, re, shutil, sys
from collections import defaultdict
from pathlib import Path

NA="NA"

def safe(s):
    s=re.sub(r"[^A-Za-z0-9_.-]+","_",str(s).strip()).strip("_")
    return s or "UNKNOWN"

def tsv_read(p):
    with open(p,newline="") as f: return list(csv.DictReader(f,delimiter="\t"))

def tsv_write(p, rows, fields):
    Path(p).parent.mkdir(parents=True,exist_ok=True)
    with open(p,"w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter="\t",extrasaction="ignore",lineterminator="\n"); w.writeheader(); w.writerows(rows)

def pair_key(name):
    return re.sub(r"_R[12](?:_001)?(?=\.f(?:ast)?q\.gz$)","",name,flags=re.I)

def parse_mlw_bio(name):
    stem=re.sub(r"\.f(?:ast)?q\.gz$","",name,flags=re.I)
    stem=re.sub(r"_R[12](?:_001)?$","",stem,flags=re.I)
    mnum=re.match(r"^(\d+)-(.+)$",stem)
    num=mnum.group(1) if mnum else NA
    x=mnum.group(2) if mnum else stem
    x=re.sub(r"_[ACGTN]+(?:-[ACGTN]+)?_L\d{3}$","",x,flags=re.I)
    x=re.sub(r"_L\d{3}$","",x,flags=re.I)
    return safe(x.replace("-","_")), num

def build_manifest(a):
    project=Path(a.project).expanduser(); raw=project/'01_raw/MLW_standard'; mlw=project/'02_analysis_fastq/MLW'; cgr=project/'02_analysis_fastq/CGR_deep'; mlw.mkdir(parents=True,exist_ok=True)
    rows=[]; pairs={}
    for f in sorted(raw.glob('*.f*q.gz')):
        key=pair_key(f.name); read='R1' if re.search(r'_R1(?:_001)?\.f',f.name,re.I) else ('R2' if re.search(r'_R2(?:_001)?\.f',f.name,re.I) else None)
        if read: pairs.setdefault(key,{})[read]=f
    prelim=[]
    for key,d in pairs.items():
        if set(d)!={'R1','R2'}: raise SystemExit(f"Unpaired MLW library: {key}: {sorted(d)}")
        bio,num=parse_mlw_bio(d['R1'].name); prelim.append((bio,num,d))
    counts=defaultdict(int)
    for bio,_,_ in prelim: counts[bio]+=1
    for bio,num,d in prelim:
        lib=f"{bio}_MLW" if counts[bio]==1 else f"{bio}_{num}_MLW"
        for read in ('R1','R2'):
            link=mlw/f"{lib}_{read}.fastq.gz"
            if link.exists() or link.is_symlink():
                if link.resolve()!=d[read].resolve(): raise SystemExit(f"Symlink collision: {link}")
            else: link.symlink_to(d[read].resolve())
        rows.append(dict(library_id=lib,biological_isolate_id=bio,sequencing_stream='MLW',original_sample_number=num,r1=str((mlw/f'{lib}_R1.fastq.gz').resolve()),r2=str((mlw/f'{lib}_R2.fastq.gz').resolve()),original_r1=str(d['R1'].resolve()),original_r2=str(d['R2'].resolve()),technical_replicate=''))
    # CGR derived FASTQs: prefer *_CGRdeep_R1; also accept *_E_R1 and normalize _E away.
    cpairs={}
    for f in sorted(cgr.glob('*.f*q.gz')):
        if '/normalised/' in str(f): continue
        read='R1' if re.search(r'_R1(?:_001)?\.f',f.name,re.I) else ('R2' if re.search(r'_R2(?:_001)?\.f',f.name,re.I) else None)
        if not read: continue
        key=pair_key(f.name); cpairs.setdefault(key,{})[read]=f
    for key,d in cpairs.items():
        if set(d)!={'R1','R2'}: raise SystemExit(f"Unpaired CGR library: {key}")
        x=key; x=re.sub(r'_CGRdeep$','',x,flags=re.I); x=re.sub(r'_E$','',x,flags=re.I); bio=safe(x.replace('-','_')); lib=f"{bio}_CGRdeep"
        # If original files are not clean-named, make clean symlinks without touching originals.
        outs={}
        for read in ('R1','R2'):
            desired=cgr/f"{lib}_{read}.fastq.gz"
            if desired != d[read]:
                if not desired.exists(): desired.symlink_to(d[read].resolve())
                outs[read]=desired
            else: outs[read]=d[read]
        rows.append(dict(library_id=lib,biological_isolate_id=bio,sequencing_stream='CGRdeep',original_sample_number=NA,r1=str(outs['R1'].resolve()),r2=str(outs['R2'].resolve()),original_r1=str(d['R1'].resolve()),original_r2=str(d['R2'].resolve()),technical_replicate=''))
    bybio=defaultdict(list)
    for r in rows: bybio[r['biological_isolate_id']].append(r)
    for r in rows: r['technical_replicate']='yes' if len(bybio[r['biological_isolate_id']])>1 else 'no'
    rows=sorted(rows,key=lambda x:x['library_id'])
    tsv_write(a.output,rows,['library_id','biological_isolate_id','sequencing_stream','original_sample_number','r1','r2','original_r1','original_r2','technical_replicate'])
    print(f"manifest libraries={len(rows)} MLW={sum(r['sequencing_stream']=='MLW' for r in rows)} CGRdeep={sum(r['sequencing_stream']=='CGRdeep' for r in rows)}")

def host_summary(a):
    man={r['library_id']:r for r in tsv_read(a.manifest)}
    rows=[]
    for lib,r in sorted(man.items()):
        mp=Path(a.host)/lib/'host_depletion_metrics.tsv'
        if not mp.exists():
            raise SystemExit(f"Missing host-depletion metrics for {lib}: {mp}")
        rr=tsv_read(mp)
        if len(rr)!=1:
            raise SystemExit(f"Expected one host-depletion metrics row for {lib}")
        m=rr[0]
        dep=float(m['retained_nonhuman_bases'])/float(a.genome_size)
        rows.append({**r,**m,'planning_depth_x':f"{dep:.3f}",
                     'normalisation_action':'downsample' if dep>float(a.cap) else 'unchanged'})
    tsv_write(a.output,rows,list(rows[0].keys()))

def fastp_summary(a):
    man=tsv_read(a.manifest); out=[]
    for r in man:
        p=Path(a.fastp)/r['library_id']/f"{r['library_id']}.json"
        with open(p) as f: j=json.load(f)
        b=j['summary']['after_filtering']['total_bases']; reads=j['summary']['after_filtering']['total_reads']
        dep=b/float(a.genome_size)
        out.append({**r,'clean_total_bases':b,'clean_total_reads':reads,'planning_depth_x':f"{dep:.3f}",'normalisation_action':'downsample' if dep>float(a.cap) else 'unchanged'})
    tsv_write(a.output,out,list(out[0].keys()))

def kraken_summary(a):
    man={r['library_id']:r for r in tsv_read(a.manifest)}
    rows=[]

    for lib,r in sorted(man.items()):
        report=Path(a.kraken)/lib/f"{lib}.report.tsv"

        if not report.is_file():
            raise SystemExit(f"Missing Kraken2 report for {lib}: {report}")

        best_species=(0,'S','NA','NA')
        best_genus=(0,'G','NA','NA')

        total_reads=None
        salmonella_clade_reads=0
        human_clade_reads=0

        parsed=[]

        with open(report) as f:
            for line in f:
                p=line.rstrip('\n').split('\t')
                if len(p)<6:
                    continue

                try:
                    pct=float(p[0])
                    clade_reads=int(p[1])
                except Exception:
                    raise SystemExit(
                        f"Unparseable Kraken2 report row for {lib}: {line.rstrip()}"
                    )

                rank=p[3].strip()
                taxid=p[4].strip()
                name=p[5].strip()

                parsed.append(
                    (pct,clade_reads,rank,taxid,name)
                )

                if rank=='U' and taxid=='0':
                    unclassified_reads=clade_reads

                if taxid=='1' and rank=='R':
                    classified_reads=clade_reads

                if taxid==str(a.salmonella_taxid):
                    salmonella_clade_reads=clade_reads

                if taxid==str(a.human_taxid):
                    human_clade_reads=clade_reads

                if rank=='S' and pct>best_species[0]:
                    best_species=(pct,rank,taxid,name)

                if rank=='G' and pct>best_genus[0]:
                    best_genus=(pct,rank,taxid,name)

        try:
            total_reads=classified_reads+unclassified_reads
        except NameError:
            raise SystemExit(
                f"Kraken2 report for {lib} lacks root and/or unclassified counts"
            )

        if total_reads<=0:
            raise SystemExit(
                f"Kraken2 report for {lib} has non-positive total read count"
            )

        sal_pct=100.0*salmonella_clade_reads/total_reads
        human_pct=100.0*human_clade_reads/total_reads

        best=best_species if best_species[0]>=10 else best_genus

        if sal_pct>=float(a.threshold):
            cls='SALMONELLA_GATE_PASS'
        elif max(best_species[0],best_genus[0])>=70:
            cls='NON_SALMONELLA_DOMINANT'
        else:
            cls='MIXED_OR_AMBIGUOUS'

        rows.append({
            **r,
            'dominant_taxon':best[3],
            'dominant_taxid':best[2],
            'dominant_rank':best[1],
            'dominant_pct':f"{best[0]:.3f}",
            'kraken_total_sequences':total_reads,
            'salmonella_clade_reads_kraken2':salmonella_clade_reads,
            'salmonella_read_pct':f"{sal_pct:.6f}",
            'residual_human_clade_reads_kraken2':human_clade_reads,
            'residual_human_read_pct_kraken2':f"{human_pct:.6f}",
            'residual_human_review_flag':'yes' if human_clade_reads>0 else 'no',
            'taxonomy_class':cls,
            'salmonella_gate_pass':'yes' if sal_pct>=float(a.threshold) else 'no',
        })

    if not rows:
        raise SystemExit("Kraken2 summary produced zero rows")

    tsv_write(a.output,rows,list(rows[0].keys()))


def parse_quast(p):
    p=Path(p)
    if not p.is_file():
        raise SystemExit(f"Missing QUAST report: {p}")

    d={}
    with open(p) as f:
        rd=csv.reader(f,delimiter="\t")
        for row in rd:
            if len(row)>=2:
                d[row[0]]=row[1]

    def require(prefix):
        for k,v in d.items():
            if k==prefix:
                try:
                    return int(float(v))
                except Exception:
                    raise SystemExit(
                        f"Unparseable QUAST value in {p}: {k}={v}"
                    )
        raise SystemExit(
            f"Required QUAST metric '{prefix}' missing from {p}"
        )

    return (
        require('# contigs'),
        require('N50'),
        require('Total length'),
    )


def qc_merge(a):
    tax={r['library_id']:r for r in tsv_read(a.taxonomy)}
    check={}

    checkm_path=Path(a.checkm)
    if not checkm_path.is_file():
        raise SystemExit(f"Missing CheckM QA table: {checkm_path}")

    with open(checkm_path,newline='') as f:
        rd=csv.DictReader(f,delimiter='\t')
        for r in rd:
            key=(r.get('Bin Id') or r.get('Bin Id ') or '').strip()
            if key:
                check[key]=r

    rows=[]

    for lib,r in sorted(tax.items()):
        q=Path(a.quast)/lib/'report.tsv'
        contigs,n50,length=parse_quast(q)

        if lib not in check:
            raise SystemExit(
                f"CheckM produced no QA row for library {lib}"
            )

        c=check[lib]

        try:
            comp=float(c['Completeness'])
            contam=float(c['Contamination'])
        except Exception:
            raise SystemExit(
                f"Could not parse CheckM completeness/contamination for {lib}"
            )

        if not math.isfinite(comp) or not math.isfinite(contam):
            raise SystemExit(
                f"Non-finite CheckM completeness/contamination for {lib}: "
                f"completeness={comp}, contamination={contam}"
            )

        reasons=[]

        if r['salmonella_gate_pass']!='yes':
            reasons.append(
                'SALMONELLA_SPECIFIC_ANALYSIS_NOT_ELIGIBLE_LT_70_PERCENT'
            )
        else:
            if contam>float(a.max_contam):
                reasons.append('CONTAMINATION_GT_20_PERCENT')
            if comp<float(a.min_comp):
                reasons.append('COMPLETENESS_LT_90_PERCENT')
            if contigs>int(a.max_contigs):
                reasons.append('CONTIGS_GT_500')
            if n50<int(a.min_n50):
                reasons.append('N50_LT_20000_BP')
            if length<int(a.min_len):
                reasons.append('GENOME_LENGTH_LT_4000000_BP')
            if length>int(a.max_len):
                reasons.append('GENOME_LENGTH_GT_5800000_BP')

        rigby=(
            'yes'
            if r['salmonella_gate_pass']=='yes' and not reasons
            else 'no'
        )

        if rigby=='yes':
            disposition='SALMONELLA_RIGBY_QC_PASS'
        elif r['salmonella_gate_pass']=='yes':
            disposition='SALMONELLA_RIGBY_QC_FAIL'
        else:
            disposition=r['taxonomy_class']

        rows.append({
            **r,
            'assembly_length_bp':length,
            'contig_count':contigs,
            'n50_bp':n50,
            'checkm_completeness':f"{comp:.3f}",
            'checkm_contamination':f"{contam:.3f}",
            'rigby_qc_pass':rigby,
            'exclusion_reasons':';'.join(reasons) if reasons else 'NONE',
            'qc_disposition':disposition,
        })

    if not rows:
        raise SystemExit("QC merge produced zero rows")

    tsv_write(a.output,rows,list(rows[0].keys()))


def typing_merge(a):
    qc={r['library_id']:r for r in tsv_read(a.qc)}; rows=[]
    for lib,r in sorted(qc.items()):
        sistr={}; sp=Path(a.typing)/lib/'sistr.tsv'
        if sp.exists():
            rr=tsv_read(sp); sistr=rr[0] if rr else {}
        mlst_scheme=mlst_st=NA; mp=Path(a.typing)/lib/'mlst.tsv'
        if mp.exists():
            line=mp.read_text().strip().split('\t');
            if len(line)>=3: mlst_scheme,mlst_st=line[1],line[2]
        genes=[]; rigbygenes=[]; plus_elements=[]; ap=Path(a.typing)/lib/'amrfinder.tsv'
        if ap.exists():
            with open(ap,newline='') as f:
                rd=csv.DictReader(f,delimiter='\t')
                for x in rd:
                    gene=x.get('Gene symbol') or x.get('Element symbol') or x.get('Name') or ''
                    element_type=(x.get('Type') or '').strip().upper()

                    if gene:
                        plus_elements.append(
                            f"{element_type}:{gene}" if element_type else gene
                        )

                    if gene and element_type=='AMR':
                        genes.append(gene)
                        ident=next((x[k] for k in x if 'Identity' in k and x[k]),None)
                        cov=next((x[k] for k in x if 'Coverage' in k and x[k]),None)

                        try:
                            if float(ident)>95 and float(cov)>95:
                                rigbygenes.append(gene)
                        except (TypeError,ValueError):
                            pass
        out={**r,
             'sistr_serovar':sistr.get('serovar',sistr.get('Serovar',NA)),
             'sistr_serogroup':sistr.get('serogroup',sistr.get('Serogroup',NA)),
             'sistr_subspecies':sistr.get('cgmlst_subspecies',sistr.get('subspecies',NA)),
             'sistr_cgmlst_st':sistr.get('cgmlst_ST',sistr.get('cgmlst_st',NA)),
             'sistr_qc_status':sistr.get('qc_status',NA),
             'sistr_qc_messages':sistr.get('qc_messages',NA),
             'mlst_scheme':mlst_scheme,'mlst_st':mlst_st,
             'amr_genes_default':','.join(sorted(set(genes))) if genes else 'NONE',
             'amr_genes_rigby_gt95_identity_coverage':','.join(sorted(set(rigbygenes))) if rigbygenes else 'NONE',
             'amrfinder_plus_elements':','.join(sorted(set(plus_elements))) if plus_elements else 'NONE'}
        rows.append(out)
    tsv_write(a.output,rows,list(rows[0].keys()))

def read_matrix(p):
    with open(p) as f:
        rd=csv.reader(f,delimiter='\t'); head=next(rd); return {row[0]:{head[i]:row[i] for i in range(1,len(head))} for row in rd}

def reconcile(a):
    rows=tsv_read(a.summary)
    by=defaultdict(list)

    for r in rows:
        by[r['biological_isolate_id']].append(r)

    mat=read_matrix(a.dist) if Path(a.dist).exists() else {}

    overrides={}
    if a.overrides and Path(a.overrides).exists():
        for x in tsv_read(a.overrides):
            overrides[x['biological_isolate_id']]=x

    rec=[]
    canon=[]

    def score(r):
        return (
            float(r['checkm_contamination']),
            -float(r['checkm_completeness']),
            int(r['contig_count']),
            -int(r['n50_bp']),
            r['library_id'],
        )

    def dval(a1,a2):
        try:
            value=mat[a1][a2]
            if value in ('','NA'):
                return None
            return float(value)
        except Exception:
            return None

    for bio,ls in sorted(by.items()):
        streams=defaultdict(list)

        for x in ls:
            streams[x['sequencing_stream']].append(x)

        ids=[x['library_id'] for x in ls]
        rigby_pass=[x for x in ls if x['rigby_qc_pass']=='yes']

        pairds=[]
        for i in range(len(ids)):
            for j in range(i+1,len(ids)):
                d=dval(ids[i],ids[j])
                if d is not None:
                    pairds.append((ids[i],ids[j],d))

        maxd=max((x[2] for x in pairds),default=None)

        same_tax=len({
            (
                x.get('taxonomy_class',NA),
                x.get('dominant_taxid',NA),
            )
            for x in ls
        })==1

        same_ser=len({
            x.get('sistr_serovar',NA)
            for x in ls
        })==1

        same_st=len({
            x.get('mlst_st',NA)
            for x in ls
        })==1

        same_amr=len({
            x.get('amr_genes_default',NA)
            for x in ls
        })==1

        salmonella_gate_count=sum(
            x.get('salmonella_gate_pass')=='yes'
            for x in ls
        )

        chosen=None
        status=None
        reason=None

        if bio in overrides:
            cid=overrides[bio]['canonical_library_id']
            chosen=next(
                (x for x in ls if x['library_id']==cid),
                None,
            )

            if not chosen:
                raise SystemExit(
                    f"Invalid replicate override for {bio}: {cid}"
                )

            status='MANUAL_OVERRIDE'
            reason=overrides[bio].get(
                'adjudication_reason',
                'manual override',
            )

        elif (
            len(streams.get('MLW',[]))>1
            or len(streams.get('CGRdeep',[]))>1
        ):
            status='REVIEW_REQUIRED'
            reason=(
                'more than one library carries the same biological label '
                'within a sequencing stream; mapping requires review'
            )

        elif len(ls)==1:
            chosen=ls[0]
            status='SINGLE_LIBRARY'
            reason='only sequencing library for this biological isolate'

        elif (
            len(ls)==2
            and len(streams.get('MLW',[]))==1
            and len(streams.get('CGRdeep',[]))==1
        ):
            if not same_tax:
                status='REVIEW_REQUIRED'
                reason=(
                    'MLW/CGRdeep representations are discordant for '
                    'dominant taxonomy; do not collapse automatically'
                )

            elif salmonella_gate_count==2:
                typing_concordant=(
                    same_ser
                    and same_st
                    and same_amr
                )

                if not typing_concordant:
                    status='REVIEW_REQUIRED'
                    reason=(
                        'both cross-stream representations pass the '
                        'Salmonella gate but are discordant for serovar, '
                        'ST, or AMR; do not collapse automatically'
                    )

                elif len(rigby_pass)==2:
                    if (
                        maxd is not None
                        and maxd<=float(a.max_tech_snp)
                    ):
                        chosen=sorted(ls,key=score)[0]
                        status=(
                            'AUTO_CONCORDANT_CROSS_STREAM_'
                            'TECHNICAL_REPLICATE'
                        )
                        reason=(
                            f'one MLW and one CGRdeep library share the '
                            f'biological label, are concordant for '
                            f'taxonomy/serovar/ST/AMR, and are <= '
                            f'{float(a.max_tech_snp):g} preliminary SNPs '
                            f'apart; canonical chosen by objective QC '
                            f'ranking. This threshold is used only for '
                            f'known cross-stream technical-replicate '
                            f'reconciliation, not transmission inference'
                        )
                    else:
                        status='REVIEW_REQUIRED'
                        reason=(
                            'both cross-stream representations pass Rigby '
                            'QC but their common-core SNP distance is '
                            'unavailable or exceeds the prespecified '
                            'technical-replicate threshold'
                        )

                elif len(rigby_pass)==1:
                    chosen=rigby_pass[0]
                    status=(
                        'CROSS_STREAM_CONCORDANT_'
                        'ONE_RIGBY_QC_PASS'
                    )
                    reason=(
                        'both cross-stream representations pass the '
                        'Salmonella gate and are concordant for '
                        'taxonomy/serovar/ST/AMR; only one passes Rigby '
                        'QC, so that representation is selected'
                    )

                else:
                    chosen=sorted(ls,key=score)[0]
                    status=(
                        'CROSS_STREAM_CONCORDANT_SALMONELLA_'
                        'NO_RIGBY_QC_PASS'
                    )
                    reason=(
                        'both cross-stream representations pass the '
                        'Salmonella gate and are concordant for '
                        'taxonomy/serovar/ST/AMR, but neither passes '
                        'Rigby QC; a QC-ranked general canonical '
                        'representation is retained but is not eligible '
                        'for the broad Salmonella tree'
                    )

            elif salmonella_gate_count==0:
                chosen=sorted(ls,key=score)[0]
                status='CROSS_STREAM_CONCORDANT_NON_SALMONELLA'
                reason=(
                    'both cross-stream representations are concordant '
                    'for dominant non-Salmonella taxonomy; canonical '
                    'representation chosen by objective general-QC '
                    'ranking'
                )

            else:
                status='REVIEW_REQUIRED'
                reason=(
                    'only one cross-stream representation passes the '
                    'Salmonella gate; do not collapse automatically'
                )

        else:
            status='REVIEW_REQUIRED'
            reason=(
                'unexpected technical-replicate configuration; '
                'manual review required'
            )

        rec.append({
            'biological_isolate_id':bio,
            'MLW_libraries':
                ';'.join(
                    x['library_id']
                    for x in streams.get('MLW',[])
                ) or NA,
            'CGRdeep_libraries':
                ';'.join(
                    x['library_id']
                    for x in streams.get('CGRdeep',[])
                ) or NA,
            'all_qc_pass_libraries':
                ';'.join(
                    x['library_id']
                    for x in rigby_pass
                ) or NA,
            'pairwise_preliminary_core_snp_distances':
                ';'.join(
                    f'{x}:{y}={d:g}'
                    for x,y,d in pairds
                ) or NA,
            'maximum_preliminary_core_snp_distance':
                NA if maxd is None else f'{maxd:g}',
            'taxonomy_concordant':'yes' if same_tax else 'no',
            'serovar_concordant':'yes' if same_ser else 'no',
            'st_concordant':'yes' if same_st else 'no',
            'amr_concordant':'yes' if same_amr else 'no',
            'status':status,
            'selected_canonical_library':
                chosen['library_id'] if chosen else NA,
            'reason':reason,
        })

        if chosen:
            z=dict(chosen)
            z['canonical_reason']=reason
            canon.append(z)

    fields=[
        'biological_isolate_id',
        'MLW_libraries',
        'CGRdeep_libraries',
        'all_qc_pass_libraries',
        'pairwise_preliminary_core_snp_distances',
        'maximum_preliminary_core_snp_distance',
        'taxonomy_concordant',
        'serovar_concordant',
        'st_concordant',
        'amr_concordant',
        'status',
        'selected_canonical_library',
        'reason',
    ]

    tsv_write(a.reconciliation,rec,fields)
    tsv_write(
        a.canonical,
        canon,
        list(rows[0].keys())+['canonical_reason'],
    )

    print(
        f"canonical={len(canon)} "
        f"review_required="
        f"{sum(x['status']=='REVIEW_REQUIRED' for x in rec)}"
    )


def candidates(a):
    metadata_path=Path(a.metadata)

    if not metadata_path.is_file():
        raise SystemExit(
            f"Required frozen isolate metadata missing: {metadata_path}"
        )

    can=[r for r in tsv_read(a.canonical) if r.get('rigby_qc_pass')=='yes']
    metadata_rows=tsv_read(metadata_path)
    meta={r['biological_isolate_id']:r for r in metadata_rows}

    missing_meta=[
        r['biological_isolate_id']
        for r in can
        if r['biological_isolate_id'] not in meta
    ]

    if missing_meta:
        raise SystemExit(
            "Canonical isolates missing from frozen isolate metadata: "
            + ",".join(sorted(set(missing_meta)))
        )

    x=[]

    for r in can:
        m=meta[r['biological_isolate_id']]
        x.append({**r,**m})

    groups=[]
    gid=0

    def add(kind,key,members):
        nonlocal gid

        if len(members)<2:
            return

        gid+=1

        for m in members:
            groups.append({
                'candidate_group_id':f'G{gid:03d}',
                'candidate_type':kind,
                'group_key':key,
                'biological_isolate_id':m['biological_isolate_id'],
                'library_id':m['library_id'],
                'serovar':m.get('sistr_serovar',NA),
                'st':m.get('mlst_st',NA),
                'household_id':m.get('household_id',NA),
                'participant_id':m.get('participant_id',NA),
                'collection_date':m.get('collection_date',NA),
                'source_type':m.get('source_type',NA),
            })

    # Same participant / same serovar / same ST.
    d=defaultdict(list)
    for r in x:
        participant=r.get('participant_id',NA)
        if participant not in ('',NA):
            d[
                (
                    participant,
                    r.get('sistr_serovar',NA),
                    r.get('mlst_st',NA),
                )
            ].append(r)

    for k,v in d.items():
        add('WITHIN_PERSON','|'.join(k),v)

    # Same household / same serovar / same ST.
    d=defaultdict(list)
    for r in x:
        household=r.get('household_id',NA)
        if household not in ('',NA):
            d[
                (
                    household,
                    r.get('sistr_serovar',NA),
                    r.get('mlst_st',NA),
                )
            ].append(r)

    for k,v in d.items():
        add('WITHIN_HOUSEHOLD','|'.join(k),v)

    # Human/environment in same household, same serovar/ST.
    d=defaultdict(list)

    for r in x:
        household=r.get('household_id',NA)
        if household not in ('',NA):
            d[
                (
                    household,
                    r.get('sistr_serovar',NA),
                    r.get('mlst_st',NA),
                )
            ].append(r)

    for k,v in d.items():
        sources={
            str(q.get('source_type','')).strip().lower()
            for q in v
        }

        if 'human' in sources and 'environment' in sources:
            add('HUMAN_ENVIRONMENT','|'.join(k),v)

    # Clinically important lineage: ST313.
    st313=[
        r for r in x
        if str(r.get('mlst_st',NA)).strip()=='313'
    ]

    add('CLINICALLY_IMPORTANT_ST313','ST313',st313)

    fields=[
        'candidate_group_id',
        'candidate_type',
        'group_key',
        'biological_isolate_id',
        'library_id',
        'serovar',
        'st',
        'household_id',
        'participant_id',
        'collection_date',
        'source_type',
    ]

    tsv_write(a.output,groups,fields)


def exclusion_ledger(a):
    rows=tsv_read(a.summary)

    rec_by_bio={}
    if Path(a.reconciliation).exists():
        rec_by_bio={
            r['biological_isolate_id']:r
            for r in tsv_read(a.reconciliation)
        }

    out=[]

    for r in rows:
        bio=r['biological_isolate_id']
        reason=r['exclusion_reasons']
        rr=rec_by_bio.get(bio)

        if rr is None:
            final='PENDING_RECONCILIATION'

        elif rr['status']=='REVIEW_REQUIRED':
            final='REPLICATE_RECONCILIATION_REVIEW_REQUIRED'
            if reason=='NONE':
                reason=rr.get('reason','REVIEW_REQUIRED')

        elif rr.get('selected_canonical_library')==r['library_id']:
            final='CANONICAL_FINAL'

        elif rr.get('selected_canonical_library') not in ('',NA):
            final='TECHNICAL_REPLICATE_NONCANONICAL'
            if reason=='NONE':
                reason='NOT_SELECTED_AFTER_TECHNICAL_REPLICATE_RECONCILIATION'

        else:
            final='UNRESOLVED_RECONCILIATION'

        out.append({
            'library_id':r['library_id'],
            'biological_isolate_id':bio,
            'sequencing_stream':r['sequencing_stream'],
            'taxonomy_class':r['taxonomy_class'],
            'salmonella_read_pct':r['salmonella_read_pct'],
            'rigby_qc_pass':r['rigby_qc_pass'],
            'qc_disposition':r['qc_disposition'],
            'exclusion_reasons':reason,
            'final_disposition':final,
        })

    if not out:
        raise SystemExit("Exclusion ledger produced zero rows")

    tsv_write(a.output,out,list(out[0].keys()))


def plot_tree(a):
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle
    from Bio import Phylo

    metadata_path=Path(a.metadata)

    if not metadata_path.is_file():
        raise SystemExit(
            f"Required frozen isolate metadata missing: {metadata_path}"
        )

    can=[
        r for r in tsv_read(a.canonical)
        if r.get('rigby_qc_pass')=='yes'
    ]
    bylib={r['library_id']:r for r in can}
    meta={
        r['biological_isolate_id']:r
        for r in tsv_read(metadata_path)
    }

    missing_meta=[
        r['biological_isolate_id']
        for r in can
        if r['biological_isolate_id'] not in meta
    ]

    if missing_meta:
        raise SystemExit(
            "Canonical isolates missing from frozen isolate metadata: "
            + ",".join(sorted(set(missing_meta)))
        )

    out=[]

    for r in can:
        m=meta[r['biological_isolate_id']]
        out.append({
            'library_id':r['library_id'],
            'biological_isolate_id':r['biological_isolate_id'],
            'subspecies':r.get('sistr_subspecies',NA),
            'serovar':r.get('sistr_serovar',NA),
            'ST':r.get('mlst_st',NA),
            'source_type':m.get('source_type',NA),
            'household_id':m.get('household_id',NA),
            'participant_id':m.get('participant_id',NA),
            'collection_date':m.get('collection_date',NA),
            'AMR':r.get(
                'amr_genes_rigby_gt95_identity_coverage',
                'NONE'
            ),
        })

    tsv_write(
        a.tree_metadata,
        out,
        [
            'library_id',
            'biological_isolate_id',
            'subspecies',
            'serovar',
            'ST',
            'source_type',
            'household_id',
            'participant_id',
            'collection_date',
            'AMR',
        ],
    )

    tree=Phylo.read(a.tree,'newick')
    n=max(1,len(tree.get_terminals()))
    fig,ax=plt.subplots(figsize=(14,max(8,n*0.22)))

    Phylo.draw(
        tree,
        axes=ax,
        do_show=False,
        show_confidence=False,
        label_func=lambda c:c.name if c.is_terminal() else None,
    )

    texts={
        t.get_text().strip():t
        for t in ax.texts
        if t.get_text().strip() in bylib
    }

    xmax=ax.get_xlim()[1]
    span=max(xmax-ax.get_xlim()[0],1e-9)
    start=xmax+0.03*span
    colw=0.035*span

    tracks=[
        ('serovar','Serovar'),
        ('ST','ST'),
        ('source_type','Source'),
        ('AMR','AMR >95/95'),
    ]

    cmap=plt.get_cmap('tab20')

    for ci,(field,label) in enumerate(tracks):
        vals=sorted({
            next(
                (
                    r[field]
                    for r in out
                    if r['library_id']==lib
                ),
                NA,
            )
            for lib in texts
        })

        colors={
            v:cmap(i%20)
            for i,v in enumerate(vals)
        }

        for lib,t in texts.items():
            y=t.get_position()[1]
            rec=next(r for r in out if r['library_id']==lib)
            v=rec[field]

            ax.add_patch(
                Rectangle(
                    (start+ci*colw,y-0.35),
                    colw*0.82,
                    0.7,
                    facecolor=colors[v],
                    edgecolor='none',
                    clip_on=False,
                )
            )

        ax.text(
            start+ci*colw+colw*0.4,
            ax.get_ylim()[1]+0.5,
            label,
            rotation=90,
            ha='center',
            va='bottom',
            fontsize=8,
            clip_on=False,
        )

    ax.set_xlim(
        ax.get_xlim()[0],
        start+len(tracks)*colw+0.04*span,
    )
    ax.set_title(
        'TiNTS culture-confirmed Salmonella diversity',
        fontsize=12,
    )
    ax.set_xlabel('Core-genome substitutions per site')
    ax.set_ylabel('')

    fig.tight_layout()
    fig.savefig(a.svg,bbox_inches='tight')
    fig.savefig(a.pdf,bbox_inches='tight')
    plt.close(fig)


def main():
    p=argparse.ArgumentParser(); sp=p.add_subparsers(dest='cmd',required=True)
    q=sp.add_parser('manifest'); q.add_argument('--project',required=True); q.add_argument('--output',required=True); q.set_defaults(fn=build_manifest)
    q=sp.add_parser('host-summary'); q.add_argument('--manifest',required=True); q.add_argument('--host',required=True); q.add_argument('--genome-size',required=True); q.add_argument('--cap',required=True); q.add_argument('--output',required=True); q.set_defaults(fn=host_summary)
    q=sp.add_parser('fastp-summary'); q.add_argument('--manifest',required=True); q.add_argument('--fastp',required=True); q.add_argument('--genome-size',required=True); q.add_argument('--cap',required=True); q.add_argument('--output',required=True); q.set_defaults(fn=fastp_summary)
    q=sp.add_parser('kraken-summary'); q.add_argument('--manifest',required=True); q.add_argument('--kraken',required=True); q.add_argument('--salmonella-taxid',required=True); q.add_argument('--human-taxid',required=True); q.add_argument('--threshold',required=True); q.add_argument('--output',required=True); q.set_defaults(fn=kraken_summary)
    q=sp.add_parser('qc'); q.add_argument('--taxonomy',required=True); q.add_argument('--quast',required=True); q.add_argument('--checkm',required=True); q.add_argument('--max-contam',required=True); q.add_argument('--min-comp',required=True); q.add_argument('--max-contigs',required=True); q.add_argument('--min-n50',required=True); q.add_argument('--min-len',required=True); q.add_argument('--max-len',required=True); q.add_argument('--output',required=True); q.set_defaults(fn=qc_merge)
    q=sp.add_parser('typing'); q.add_argument('--qc',required=True); q.add_argument('--typing',required=True); q.add_argument('--output',required=True); q.set_defaults(fn=typing_merge)
    q=sp.add_parser('reconcile'); q.add_argument('--summary',required=True); q.add_argument('--dist',required=True); q.add_argument('--max-tech-snp',required=True); q.add_argument('--overrides'); q.add_argument('--reconciliation',required=True); q.add_argument('--canonical',required=True); q.set_defaults(fn=reconcile)
    q=sp.add_parser('candidates'); q.add_argument('--canonical',required=True); q.add_argument('--metadata',required=True); q.add_argument('--output',required=True); q.set_defaults(fn=candidates)
    q=sp.add_parser('ledger'); q.add_argument('--summary',required=True); q.add_argument('--reconciliation',required=True); q.add_argument('--output',required=True); q.set_defaults(fn=exclusion_ledger)
    q=sp.add_parser('plot-tree'); q.add_argument('--tree',required=True); q.add_argument('--canonical',required=True); q.add_argument('--metadata',required=True); q.add_argument('--tree-metadata',required=True); q.add_argument('--svg',required=True); q.add_argument('--pdf',required=True); q.set_defaults(fn=plot_tree)
    a=p.parse_args(); a.fn(a)
if __name__=='__main__': main()
