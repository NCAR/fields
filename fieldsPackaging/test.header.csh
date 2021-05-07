
foreach i ( *.Rd )
#echo $i
    sed  "1,3 d" $i > hold
    head -1 hold
end


