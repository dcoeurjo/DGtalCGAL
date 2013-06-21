#!/bin/zsh

typeset file=$1
typeset myH=1.0
typeset minH=0.1
typeset myR=6.0



while (( myH> minH ))
do
    echo " =========== h $myH ==========="

    echo "Curvatures computation..."
    3dLocalEstimators "x^2+y^2+z^2-25" -10 -10 -10 10 10 10  $myH $myR >! CurvatureCircle_$myH.txt


    echo "Computing Stats..."
    statisticsEstimators CurvatureCircle_$myH.txt 5 9 | tail -1 >> ErrorMongeMean.txt
    Statisticsestimators CurvatureCircle_$myH.txt 6 10 | tail -1 >> ErrorMongeGaussian.txt
    statisticsEstimators CurvatureCircle_$myH.txt 5 7 | tail -1 >> ErrorIIMean.txt
    statisticsEstimators CurvatureCircle_$myH.txt 6 8 | tail -1 >> ErrorIIGaussian.txt

    ### Crappy but it's just to get first results
    let "myH=$myH - 0.1"

done    
