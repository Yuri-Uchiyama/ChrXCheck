all:
	javac jp/ac/yokohama_cu/genetics/ngs/ChrXCheck.java
test:
	java jp.ac.yokohama_cu.genetics.ngs.ChrXCheck -s depth.00.sample_interval_summary -n normalized.sample_interval_summary -r list.txt -d deduced.txt
help:
	java jp.ac.yokohama_cu.genetics.ngs.ChrXCheck 
