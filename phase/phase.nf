






/mnt/ebs/ref_push/capture/HLA-1/bam/
PC-HLA-EG4.bam      PC-HLA-EG5.bam      PC-HLA-EG6.bam      PC-HLA-GM-1.bam      PC-HLA-GM-2.bam      PC-HLA-GM-4.bam
PC-HLA-EG4.bam.bai  PC-HLA-EG5.bam.bai  PC-HLA-EG6.bam.bai  PC-HLA-GM-1.bam.bai  PC-HLA-GM-2.bam.bai  PC-HLA-GM-4.bam.bai

/mnt/ebs/james/hla/capture/dvout/
PC-HLA-EG4.g.vcf.gz  PC-HLA-EG5.g.vcf.gz  PC-HLA-EG6.g.vcf.gz  PC-HLA-GM-1.g.vcf.gz  PC-HLA-GM-2.g.vcf.gz  PC-HLA-GM-4.g.vcf.gz
PC-HLA-EG4.vcf.gz    PC-HLA-EG5.vcf.gz    PC-HLA-EG6.vcf.gz    PC-HLA-GM-1.vcf.gz    PC-HLA-GM-2.vcf.gz    PC-HLA-GM-4.vcf.gz

extractHAIRS --indels 1 --hic 1 --bam c6/c6_40X_HG002.bam --VCF c6/c6_40X_variants_clean.vcf --out c6_withindels/c6_40X_HG002.fragments
HAPCUT2 --hic 1 --fragments c6_withindels/c6_40X_HG002.fragments --VCF c6/c6_40X_variants_clean.vcf --output c6_withindels/c6_40X_HG002_hapcut_output --outvcf 1 




extractHAIRS --indels 1 --hic 1 --bam bamfile --VCF vcfin --out fragments
HAPCUT2 --hic 1 --fragments fragments --VCF vcfin --output hapcut_output --outvcf 1 
