from __future__ import division
import sys

"""
Find the closest sequences to each antibody tag.
>>> python antibody_tag_search.py MPLELTQSRVQKIWVPVDHRPSLPRSCGPKLTNSPTVIVMVGLPARGKTYISKKLT
"""


tag_dict = dict(flag='DYKDDDDK', ha='YPYDVPDYA', his='HHHHHH', myc='EQKLISEEDL', v5='GKPIPNPLLGLDST',
                xpress='DLDDDDK', thrombin='LVPRGS')

def find_tag(target):
    for key, tag in tag_dict.iteritems():
        score = []
        for num, ti in enumerate(target):
            try:
                seq = target[num:num + len(tag)]
                score.append(len([i for i, ii in zip(seq, tag) if i == ii]))
            except:
                pass
        print key, max(score), max(score)/len(tag), score.index(max(score))


if __name__ == "__main__":
    find_tag(sys.argv[1])
