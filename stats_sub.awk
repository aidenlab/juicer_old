{
a1+=$1;
a2+=$2;
a3+=$3;
a4+=$4;
a5+=$5;
}
END{
	printf("Total:  %'d\n Unmapped: %'d (%0.2f%)\n Regular: %'d (%0.2f%)\n Normal chimeric: %'d (%0.2f%) \n Abnormal chimeric: %'d (%0.2f%)\n Total alignable reads: %'d (%0.2f%)\n", a1, a2, a2*100/a1, a3, a3*100/a1, a4, a4*100/a1, a5, a5*100/a1, a3+a4, (a3+a4)*100/a1);
}
