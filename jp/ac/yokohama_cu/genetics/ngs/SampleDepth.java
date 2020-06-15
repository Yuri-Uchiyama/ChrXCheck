package jp.ac.yokohama_cu.genetics.ngs;

public class SampleDepth {

	int[] auto_total;
	double[] auto_mean;
	int[] auto_gq1;
	int[] auto_med;
	int[] auto_gq3;
	double[] auto_a15;

	int[] x_total;
	double[] x_mean;
	int[] x_gq1;
	int[] x_med;
	int[] x_gq3;
	double[] x_a15;

	int[] el_total;
	double[] el_mean;
	int[] el_gq1;
	int[] el_med;
	int[] el_gq3;
	double[] el_a15;


	int autoCounter = 0;
	int xCounter = 0;
	int elseCounter = 0;

	public SampleDepth(int sizeofAutosomal, int sizeofXChromosomal, int sizeofElseChromosomal){
		auto_total = new int[sizeofAutosomal];
		auto_mean = new double[sizeofAutosomal];
		auto_gq1 = new int[sizeofAutosomal];
		auto_med = new int[sizeofAutosomal];
		auto_gq3 = new int[sizeofAutosomal];
		auto_a15 = new double[sizeofAutosomal];
		x_total = new int[sizeofXChromosomal];
		x_mean = new double[sizeofXChromosomal];
		x_gq1 = new int[sizeofXChromosomal];
		x_med = new int[sizeofXChromosomal];
		x_gq3 = new int[sizeofXChromosomal];
		x_a15 = new double[sizeofXChromosomal];

		el_total = new int[sizeofElseChromosomal];
		el_mean= new double[sizeofElseChromosomal];
		el_gq1 = new int[sizeofElseChromosomal];
		el_med = new int[sizeofElseChromosomal];
		el_gq3 = new int[sizeofElseChromosomal];
		el_a15 = new double[sizeofElseChromosomal];
	}
	public void add(String region, String total, String mean, String gq1, String med, String gq3, String a15){
		if(region.matches("^[1-9].*")){
			auto_total[autoCounter] = Integer.parseInt(total);
			auto_mean[autoCounter] = Double.parseDouble(mean);
			auto_gq1[autoCounter] = Integer.parseInt(gq1);
			auto_med[autoCounter] = Integer.parseInt(med);
			auto_gq3[autoCounter] = Integer.parseInt(gq3);
			auto_a15[autoCounter] = Double.parseDouble(a15);
			autoCounter++;
		}else if(region.startsWith("X")){
			x_total[xCounter] = Integer.parseInt(total);
			x_mean[xCounter] = Double.parseDouble(mean);
			x_gq1[xCounter] = Integer.parseInt(gq1);
			x_med[xCounter] = Integer.parseInt(med);
			x_gq3[xCounter] = Integer.parseInt(gq3);
			x_a15[xCounter] = Double.parseDouble(a15);
			xCounter++;
		}else {
			el_total[elseCounter] = Integer.parseInt(total);
			el_mean[elseCounter] = Double.parseDouble(mean);
			el_gq1[elseCounter] = Integer.parseInt(gq1);
			el_med[elseCounter] = Integer.parseInt(med);
			el_gq3[elseCounter] = Integer.parseInt(gq3);
			el_a15[elseCounter] = Double.parseDouble(a15);
			elseCounter++;
		}
	}

	// 
	// Sample_24116_total_cvg
	// Sample_24116_mean_cvg
	// Sample_24116_granular_Q1
	// Sample_24116_granular_median
	// Sample_24116_granular_Q3
	// Sample_24116_%_above_15
}
