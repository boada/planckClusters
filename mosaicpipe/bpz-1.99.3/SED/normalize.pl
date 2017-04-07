#!/usr/local/bin/perl

use PDL;

@files = glob("bpgs*.sed");


@files = ("sey1.sed",
	  "sey2.sed");

$total = 10000;
foreach $file (@files){

    ($w,$f) = rcols($file,0,1);
    
    # Fix negative values of flux to be zero
    ($tmp =  $f->where ($f < 0)) .= 0.0;    

    
    $W = $w->where($f > 0 );    
    $F = $f->where($f > 0 );    


    $dx = ($W - rotate($W,1));
    $dx = $dx->slice('1:-1');
    $y  = $F->slice('1:-1');
	
    $ftotal = intover($dx*$y);
    $scale  = $total/$ftotal;

    printf "%s %e %e\n",$file,$ftotal,$scale;

    open(OUT,">$file");
    printf OUT "%8.2f %e\n",0.00,0.00;
    wcols "%8.2f %e",$W,$F*$scale, *OUT;
    close (OUT);
    
}





