#!/bin/tcsh
cat  > cmdfile <<EOF
 s/Copyright 2004-2011/Copyright 2004-2013/g
EOF

 
foreach i ( *.R )
 echo $i
 mv $i $i.orig
 sed -f cmdfile $i.orig > $i
end

foreach i ( *.Rd )
 echo $i
 mv $i $i.orig
 sed -f cmdfile $i.orig > $i
end

