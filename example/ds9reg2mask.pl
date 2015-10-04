#!/usr/bin/perl
#
# convert ds9 region file to pixel mask file
# NOTE: support only `box' and `ellipse', ignore ang!
# NOTE: old verison of ds9 may not work
#
# usage: > ds9reg2mask.pl hoge.reg nxpix nypix
#
# output file can be converted to mask fits file using iraf:
#  > text2mask hoge.dat hoge.pl nxpix nypix
#  > imcopy hoge.pl hoge.fits
#

if(($#ARGV+1)!=3){
 die "usage:\n > ds9reg2mask.pl hoge.reg nxpix nypix\n";
}

open(IN,"$ARGV[0]") || die "File Not Found\n";
$nx=$ARGV[1];
$ny=$ARGV[2];

$nb=0;
$ne=0;
while(<IN>){
    if(/,/){
#	$n1=index($_,";");
	$n1=-1;
	$n2=index($_,"(");
	$n3=index($_,")");
	$name=substr($_,$n1+1,$n2-$n1-1);
	$para=substr($_,$n2+1,$n3-$n2-1);
	@data=split(/,/,$para);

	if($name eq "box"){
	    $nb++;
	    $bx1[$nb]=$data[0]-0.5*$data[2];
	    $bx2[$nb]=$data[0]+0.5*$data[2];
	    $by1[$nb]=$data[1]-0.5*$data[3];
	    $by2[$nb]=$data[1]+0.5*$data[3];
	}
	
	if($name eq "ellipse"){
	    $ne++;
	    $ecx[$ne]=$data[0];
	    $ecy[$ne]=$data[1];
	    $erx[$ne]=$data[2];
	    $ery[$ne]=$data[3];
	}
    }
}

close(IN);

for($i=1;$i<=$nx;$i++){
    for($j=1;$j<=$ny;$j++){
	$f=0;
	for($k=1;$k<=$nb;$k++){
	    if(($i>=$bx1[$k])&&($i<=$bx2[$k])&&($j>=$by1[$k])&&($j<=$by2[$k])){
		$f=1;
	    }
	}
	for($k=1;$k<=$ne;$k++){
	    $rr=($i-$ecx[$k])*($i-$ecx[$k])/($erx[$k]*$erx[$k])+($j-$ecy[$k])*($j-$ecy[$k])/($ery[$k]*$ery[$k]);
	    if($rr<=1.0){
		$f=1;
	    }
	}
	if($f==1){ printf "%d %d\n",$i,$j; }
    }
}


