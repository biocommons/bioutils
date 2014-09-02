import base64
import hashlib


def seq_seguid(seq):
    """returns BioPython-compatible seguid""" 
    seq = seq.upper().rstrip('*')
    return base64.standard_b64encode(hashlib.sha1(seq).digest()).rstrip('=')

def seq_md5(seq):
    "returns sequence md5 as hex digest"
    seq = seq.upper().rstrip('*')
    return hashlib.md5(seq).hexdigest()

def seq_sha1(seq):
    """returns sequence sha1 as url-safe base64 encoded digest. This is
    identical to seguid except for hashes that contain / or +."""
    seq = seq.upper().rstrip('*')
    return base64.urlsafe_b64encode(hashlib.sha1(seq).digest()).rstrip('=')



def tx_digest(seq_md5, cds_se_i, exons_se_i):
    def coord_fmt(se):
        return "[{se[0]};{se[1]})".format(se=se)
    tx_info = "{seq_md5};{cds_se_str};[{exon_se_str}]".format(
        seq_md5=seq_md5, cds_se_str=coord_fmt(cds_se_i),
        exon_se_str=";".join([coord_fmt(ex) for ex in sorted(exons_se_i)])
        )
    return tx_info
