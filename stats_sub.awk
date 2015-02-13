{
a1+=$2; # total reads
a2+=$3; # normal
a3+=$4; # chimeric paired
a4+=$5; # chimeric ambiguous
a5+=$6; # unmapped
a6+=$1; # ligations 
}
END{
    printf("Sequenced Read Pairs:  %'d\n Normal Paired: %'d (%0.2f%)\n Chimeric Paired: %'d (%0.2f%)\n Chimeric Ambiguous: %'d (%0.2f%)\n Unmapped: %'d (%0.2f%)\n Ligation Motif Present: %'d (%0.2f%)\nAlignable (Normal+Chimeric Paired): %'d (%0.2f%)\n", a1, a3, a3*100/a1, a4, a4*100/a1, a5, a5*100/a1, a2, a2*100/a1, a6, a6*100/a1, a3+a4, (a3+a4)*100/a1);
}
