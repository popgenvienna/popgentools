import sys
import random

reps=set([(1,2),(3,4),(5,6),(7,8)])
xps=set([(1,3),(2,4),(1,5),(2,6),(1,7),(2,8),(3,5),(4,6),(3,7),(4,8),(5,7),(6,8)])
freqclasses=(0.2, 0.225, 0.25, 0.275,
             0.3, 0.325, 0.35, 0.375,
             0.4, 0.425, 0.45, 0.475,
             0.5, 0.525, 0.55, 0.575,
             0.6, 0.625, 0.65, 0.675,
             1.0)


"""
(1, 2)	0.25
(1, 2)	0.25
(1, 2)	0.25
(1, 2)	0.25
(1, 2)	0.25
(1, 2)	0.25
(1, 2)	0.25
"""

ra= [0 for i in range(0,len(freqclasses))]
xa= [0 for i in range(0,len(freqclasses))]
totr=0
totx=0

for l in open(sys.argv[1]):
    l=l.rstrip("\n")
    a=l.split("\t")
    if l.startswith("@totcounts"):
        """
        @totcounts	(1, 2)	444824
        @totcounts	(7, 8)	448141
        """
        key=eval(a[1])
        if key in reps:
            totr+=int(a[2])
        elif key in xps:
            totx+=int(a[2])
        else:
            raise Exception("invalid key")

    if not l.startswith("("):
        continue
    
    val=float(a[1])
    counter=0;
    for v in freqclasses:
        if val< v:
            break
        counter+=1
    

    key=eval(a[0])
    if key in reps:
        ra[counter]+=1
    elif key in xps:
        xa[counter]+=1
    else:
        raise Error("invalid state")

for i,f in enumerate(freqclasses):
    vr=ra[i]
    print "rep\t{0}\t{1}\t{2}\t{3}".format(f,vr,totr,float(vr)/float(totr))

for i,f in enumerate(freqclasses):
    vx=xa[i]
    print "xp\t{0}\t{1}\t{2}\t{3}".format(f,vx,totx,float(vx)/float(totx))
