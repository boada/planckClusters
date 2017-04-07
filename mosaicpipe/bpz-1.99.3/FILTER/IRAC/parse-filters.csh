# Responses from
# http://irsa.ipac.caltech.edu/data/SPITZER/docs/irac/calibrationfiles/spectralresponse/
awk '{if ($1!="#") printf "%.1f   %s\n",$1*10000,$2}' 080924ch1trans_full.txt > CH1.res
awk '{if ($1!="#") printf "%.1f   %s\n",$1*10000,$2}' 080924ch2trans_full.txt > CH2.res
