BEGIN{
	tot=0;
	name=0;
	waitstring="";
}
{
	if (tot >= 1000000) {
		if (p1 != $1 || p2 != $2 || p4 != $4 || p5 != $5 || p8 != $8) {
			sname=sprintf("%s_msplit%04d_", groupname, name);
			sysstring = sprintf("bsub -o %s -q %s -g %s -J %s \"awk -f /broad/aidenlab/neva/neva_scripts/dups.awk -v name=%s/%s %s/split%d; \"", outfile, queue, groupname, sname, dir, sname, dir, name, dir, name);
			system(sysstring);
			if (name==0) {
				waitstring=sprintf("exit(%s)", sname);
			}
			else {
				waitstring=sprintf("%s || exit(%s)", waitstring, sname);
			}
			name++;
			tot=0;
		}
	}
	outname = sprintf("%s/split%d", dir, name);
	print > outname;
	p1=$1;p2=$2;p4=$4;p5=$5;p6=$6;p8=$8;
	tot++;
}
END {
	sname=sprintf("%s_msplit%04d_", groupname, name);
	sysstring = sprintf("bsub -o %s -q %s -g %s -J %s \"awk -f /broad/aidenlab/neva/neva_scripts/dups.awk -v name=%s/%s %s/split%d; \"", outfile, queue, groupname, sname, dir, sname, dir, name);
	system(sysstring);
	if (name==0) {
		waitstring=sprintf("exit(%s)", sname);
	}
	else {
		waitstring=sprintf("%s || exit(%s)", waitstring, sname);
	}
	sysstring = sprintf("bsub -o %s -q %s -g %s -J %s_split -w \"done(%s_msplit*)\" \"cat %s/*_msplit*_optdups.txt > %s/opt_dups.txt;  cat %s/*_msplit*_dups.txt > %s/dups.txt;cat %s/*_msplit*_merged_nodups.txt > %s/merged_nodups.txt; \" ", outfile, queue, groupname, groupname, groupname, dir, dir, dir, dir, dir, dir);
	system(sysstring);
	sysstring = sprintf("bsub -o %s -q %s -g %s_skill -w \"%s\" \"bkill -g %s 0; bkill -g %skill 0\" ",outfile, queue, groupname, waitstring, groupname, groupname);
	system(sysstring);
 	sysstring = sprintf("bsub -o %s -q %s -g %s -w \"done(%s_split)\" \"bkill -g %s_skill 0 ; rm %s/*_msplit*_optdups.txt; rm %s/*_msplit*_dups.txt; rm %s/*_msplit*_merged_nodups.txt; rm %s/split* \"", outfile, queue, groupname, groupname, groupname, dir, dir, dir, dir);
	system(sysstring);


}
