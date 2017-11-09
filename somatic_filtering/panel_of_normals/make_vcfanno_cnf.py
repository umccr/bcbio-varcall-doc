import os

fpaths = '''MH17B001P004-germline-ensemble-annotated.vcf.gz
MH17B001P013-germline-ensemble-annotated.vcf.gz
VPH52_Blood-germline-ensemble-annotated.vcf.gz
VPH54_Blood-germline-ensemble-annotated.vcf.gz
VPH56_Blood-germline-ensemble-annotated.vcf.gz
VPH58_Blood-germline-ensemble-annotated.vcf.gz
VPH59_Blood-germline-ensemble-annotated.vcf.gz
VPH61_Blood-germline-ensemble-annotated.vcf.gz
WPT-013-normal-ensemble-annotated.vcf.gz
'''.split()


def lbl_by_path(fp):
    fname = os.path.basename(fp)
    fname = fname.replace('-germline-ensemble-annotated.vcf.gz', '')
    fname = fname.replace('-ensemble-annotated.vcf.gz', '')
    fname = fname.replace('.vcf.gz', '')
    fname = fname.replace('-', '_')
    return 'PoN_' + fname


for fp in fpaths:
    print(f'''[[annotation]]
file="/home/vlad/validation/normals/{fp}"
fields=[""]
ops=["self"]
names=["{lbl_by_path(fp)}"]
'''
)

print(f'''[[postannotation]]
name="PoN_CNT"
fields=[{', '.join('"' + lbl_by_path(fp) + '"' for fp in fpaths)}]
op="lua:count_true({', '.join(lbl_by_path(fp) for fp in fpaths)})"
type="Integer"

[[postannotation]]
fields=[{', '.join('"' + lbl_by_path(fp) + '"' for fp in fpaths)}]
op="delete"
'''
)
