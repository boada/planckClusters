#!/usr/local/bin/perl

use PDL;
use PDL::Graphics::PGPLOT;
use PDL::Graphics::PGPLOT::Window;
use PGPLOT;

@files = glob("bpgs*.sed");

$total = 10000;
foreach $file (@files){

    ($w,$f) = rcols($file,0,1);
    
    
    ($root,$nada) = split('\.',$file);

    ($nada,$number) = split('_',$root);

    if ($number < 10){
	$number = "00$number";
    }elsif($number < 100){ 	
	$number = "0$number";
    }


    $dev = "$number.ps/ps";
    #$dev = "/xs";

    ($x1,$x2) = minmax($w);
    ($y1,$y2) = minmax($f);
    #new(Device => '/xs');
    #pgsci(1);
    $win = dev ($dev,
		{AxisColour=>1,
		 XTitle=>"Lambda",
		 YTitle=>"Flux",
		 Title=>$file,
		 COLOR=>1});
    
    print "$file $root $number \n";

    #env($x1,$x2,$y1,$y2,
    #{COLOR=>1});
    line($w,$f);
    #hold;

    # Fix negative values of flux to be zero
    #($tmp =  $f->where ($f < 0)) .= 0.0;    

    #$dx = ($w - rotate($w,1));
    #$dx = $dx->slice('1:-1');
    #$y  = $f->slice('1:-1');
	
    #$F = intover($dx*$y);
    #$scale = $total/$F;
    #printf "%s %e %e\n",$file,$F,$scale;


    #wcols "%6.2f %e",$w,$f*$scale, $file;

}





