import sys
vcf_file = sys.argv[1]


# Parse REF, ALT, AF, MSI, MSILEN
# Add MSI_FAIL


def proc_line(ref, alt, af, msi, msilen):
    # Filter low AF MSI
    if msi:
        msi = float(msi)
        change_len = abs(len(ref) - len(alt))
        if change_len == 1 and msi > 1:
            msi_fail = any([
                msi <=  2 and af < 0.005,
                msi <=  4 and af < 0.01,
                msi <=  7 and af < 0.03,
                msi ==  8 and af < 0.06,
                msi ==  9 and af < 0.125,
                msi == 10 and af < 0.175,
                msi == 11 and af < 0.25,
                msi == 12 and af < 0.3,
                msi >  12 and af < 0.35])
            if msi_fail:
                return False
        elif change_len == 3 and msi >= 5 and af < 0.1:  # ignore low AF in 3nt MSI region
            return False
    return True


def use_pyvcf():
    """ Working, but doesn't write header.
        Time:
    """
    import vcf
    import gzip
    with (gzip.open(vcf_file) if vcf_file.endswith('.gz') else open(vcf_file)) as f:
        vcf_reader = vcf.Reader(f)
        vcf_writer = vcf.Writer(sys.stdout, vcf_reader)
        for rec in vcf_reader:
            msi_fail = proc_line(rec.REF, rec.ALT[0], rec.samples[0]['AF'], rec.INFO['MSI'], rec.INFO['MSILEN'])
            rec.FILTER.append('MSI_FAIL')
            vcf_writer.write_record(rec)

def use_cyvcf():
    import vcf  # need to reinstall instead of pyvcf
    import gzip
    with (gzip.open(vcf_file) if vcf_file.endswith('.gz') else open(vcf_file)) as f:
        vcf_reader = vcf.Reader(f)
        vcf_writer = vcf.Writer(sys.stdout, vcf_reader)
        for rec in vcf_reader:
            msi_fail = proc_line(rec.REF, rec.ALT[0], rec.samples[0]['AF'], rec.INFO['MSI'], rec.INFO['MSILEN'])
            rec.FILTER.append('MSI_FAIL')
            vcf_writer.write_record(rec)

def use_pysam():
    """ Working.
        Time:
    """
    import pysam
    vcf = pysam.VariantFile(vcf_file)
    vcf.header.filters.add('MSI_FAIL', None, None, '')
    sys.stdout.write(str(vcf.header))
    for rec in vcf:
        msi_fail = proc_line(rec.ref, rec.alts[0], rec.samples.values()[0]['AF'], rec.info['MSI'], rec.info['MSILEN'])
        rec.filter.add('MSI_FAIL')
        sys.stdout.write(str(rec))

def use_cyvcf2():
    """ segmentation fault
    """
    from cyvcf2 import VCF
    for rec in VCF(vcf_file):
        msi_fail = proc_line(rec.REF, rec.ALT[0], rec.format('AF')[0], rec.INFO['MSI'], rec.INFO['MSILEN'])
        if rec.FILTER:
            rec.FILTER += ';MSI_FAIL'
        else:
            rec.FILTER = 'MSI_FAIL'
        print(rec)

def use_python():
    """ Working.
        Time:
    """
    import re
    import gzip
    with (gzip.open(vcf_file) if vcf_file.endswith('.gz') else open(vcf_file)) as f:
        for l in f:
            if l.startswith('#'):
                sys.stdout.write(l)
            else:
                c, p, i, ref, alt, q, filt, info, f = l.split('\t')[:9]
                samples = l.split('\t')[9:]
                m = re.match(r'.*;?MSI=(?P<msi>\d);MSILEN=(?P<msi_len>\d);?.*', info)
                if m:
                    msi = m.group('msi')
                    msi_len = m.group('msi_len')
                else:
                    sys.stdout.write(l)
                af = float(samples[0].split(':')[3])
                msi_fail = proc_line(ref, alt, af, msi, msi_len)
                filt = ';'.join(filt.split(';') + [';MSI_FAIL'])
                sys.stdout.write('\t'.join([c, p, i, ref, alt, q, filt, info, f] + samples))


# use_python()
# use_pysam()
use_pyvcf()


def call_somatic():
    new_headers = ['##INFO=<ID=SOMATIC,Number=0,Type=Flag,Description="Somatic event">\n',
                   ('##FILTER=<ID=REJECT,Description="Not somatic due to normal call frequency '
                    'or phred likelihoods: tumor: %s, normal %s.">\n')
                   % (int(tumor_thresh * 10), int(normal_thresh * 10))]
    def _output_filter_line(line, indexes):
        parts = line.split("\t")
        if _check_lods(parts, tumor_thresh, normal_thresh, indexes) and _check_freqs(parts, indexes):
            parts[7] = parts[7] + ";SOMATIC"
        else:
            if parts[6] in set([".", "PASS"]):
                parts[6] = "REJECT"
            else:
                parts[6] += ";REJECT"
        line = "\t".join(parts)
        sys.stdout.write(line)
    def _write_header(header):
        for hline in header[:-1] + new_headers + [header[-1]]:
            sys.stdout.write(hline)
    header = []
    indexes = None
    for line in sys.stdin:
        if not indexes:
            if line.startswith("#"):
                header.append(line)
            else:
                parts = header[-1].rstrip().split("\t")
                indexes = {"tumor": parts.index(tumor_name), "normal": parts.index(normal_name)}
                _write_header(header)
                _output_filter_line(line, indexes)
        else:
            _output_filter_line(line, indexes)
    # no calls, only output the header
    if not indexes:
        _write_header(header)

