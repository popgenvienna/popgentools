import sys
import random
from optparse import OptionParser, OptionGroup

"""
Dr. Robert Kofler
INPUT
2L	15	A	A/C	.		670:1330	527:1473	473:1527	450:1550	456:1544	556:1444	449:1551	257:1743	751:1249	820:1180	877:1123	1024:976	1084:916	1094:906	1268:732	702:1298	638:1362	663:1337	759:1241	918:1082	1047:953	1032:968	741:1259	832:1168	689:1311	731:1269	609:1391	567:1433	697:1303	871:1129	979:1021	1115:885	1223:777	1150:850	1150:850	1176:824	525:1475	495:1505	571:1429	745:1255	741:1259	646:1354	794:1206	800:1200	660:1340	692:1308	639:1361	693:1307	798:1202	721:1279	898:1102	878:1122	886:1114	860:1140	946:1054	1079:921	893:1107	590:1410	694:1306	759:1241	790:1210	662:1338	702:1298	944:1056	588:1412	680:1320	655:1345	629:1371	684:1316	457:1543	680:1320

OUTPUT
Unknown_group_104      5943    N       0:0:10:0:0:0    0:0:10:0:0:0    0:0:10:0:0:0
Unknown_group_104      5944    N       0:8:0:0:0:0     0:8:0:0:0:0     0:8:0:0:0:0

"""

def convert_entry(entry,coverage):
    a1,a2=map(float,entry.split(":"))
    sum=a1+a2
    na1=int((a1*coverage)/sum)
    na2=int(coverage)-na1
    return "{0}:0:{1}:0:0:0".format(na1,na2)

parser = OptionParser()
parser.add_option("--sum",dest="sum",help="the summary file")
parser.add_option("--coverage",dest="coverage",help="the coverage to subsample")
(options, args) = parser.parse_args()

cov=float(options.coverage)
ifh=open(options.sum)
for line in ifh:
    line=line.rstrip()
    a=line.split("\t")
    a.pop(3) # remove A/C
    a.pop(3) # remove . comment
    if(a[3]==""):
        a.pop(3) # remove bug-column..
    for i in range(3,len(a)):
        e=a[i]
        entry=convert_entry(e,cov)
        a[i]=entry
    topr="\t".join(a)
    print topr
        


