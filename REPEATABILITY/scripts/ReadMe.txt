

logic:
keep R function which gives redundant list of repeats (I allow motivs overlap each other - simple and error prone backbone).
After this there are several steps of filtration:
=> if mismatch is at the end of motive - remove this part
etc...

when I am OK with this backbone - run it for whole mtDNA (16500*3) on R server

After - I save such files for the whole mtDNA to keep the metadata - using such metadata we can fastly do statistics, visualise, etc,,,


After I use more filters to keep
only repeats from major arc, only in contact zone, only with degradation < 5% for example, only with high GC content (make some metric....)
for some of these conditions I may get nice results - it can be like machine learning - find filtration rules based on GC content, length, mismatch fraction.
I have to reread bacterial (old cell papers) to derive math and after - discuss with Falkenberg this. 

To Do Next:
=> check that redundant output is correct (do control by hand and by logic - have universal R function with robust code - minimum libraries)
=> add to redudant outputs naive number of hydrigen bonds (GC=3, AT=2), % of mismatch and total length. It should be enough to understand main effects and after - correlate with Nathan's results. Add maybe total local length (where there is a SNP)

