#!/usr/bin/env python
"""Build transcript-to-gene mapping from GTF file."""
import gzip
import re
import sys

def parse_gtf_attributes(attr_string):
    """Parse GTF attribute string into dict."""
    attrs = {}
    for attr in attr_string.strip().split(';'):
        attr = attr.strip()
        if not attr:
            continue
        match = re.match(r'(\w+)\s+"([^"]*)"', attr)
        if match:
            attrs[match.group(1)] = match.group(2)
    return attrs

def main():
    gtf_path = sys.argv[1]
    out_path = sys.argv[2]

    seen = set()

    with gzip.open(gtf_path, 'rt') as f_in, open(out_path, 'w') as f_out:
        for line in f_in:
            if line.startswith('#'):
                continue

            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue

            feature_type = fields[2]
            if feature_type != 'transcript':
                continue

            attrs = parse_gtf_attributes(fields[8])

            tid = attrs.get('transcript_id', '')
            gid = attrs.get('gene_id', '')
            gname = attrs.get('gene_name', '')

            if tid and gid and tid not in seen:
                f_out.write(f"{tid}\t{gid}\t{gname}\n")
                seen.add(tid)

    print(f"Wrote {len(seen)} transcript-gene mappings to {out_path}")

if __name__ == "__main__":
    main()
