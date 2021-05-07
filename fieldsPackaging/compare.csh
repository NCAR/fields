#!/bin/tcsh
# Comparing directories
set DIR1 = /Users/nychka/Home/Src/R.new
set DIR2 = /Users/nychka/Home/Src/R.old
set HOLD = /Users/nychka/Home/Src/LatticeKrigDiffs

cd $DIR1
 pwd 
foreach i ( *.R )
 echo $i
 diff $i $DIR2/$i >  $HOLD/$i.diff 
 diff $i $DIR2/$i | wc -l 
end

foreach i ($DIR1/$*.R )
 echo $i
 diff $i $DIR2/$i >  $HOLD/$i.diff 
 diff $i $DIR2/$i | wc -l 
end
