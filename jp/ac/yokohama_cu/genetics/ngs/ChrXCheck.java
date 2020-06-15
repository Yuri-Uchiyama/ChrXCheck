package jp.ac.yokohama_cu.genetics.ngs;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.HashMap;

public class ChrXCheck {
	public class UsedInfo {
		String name;
		String gender;
		double fvalue;
		double nfvalue;
	};
	HashMap<String, String> genders = new HashMap<>();
	HashMap<String, UsedInfo> usedGenders = new HashMap<>();
	HashMap<String, SampleDepth> depth = new HashMap<>();
	boolean normalize = false;
	boolean skipNonExistentSamples = false;
	String[] target = null;
	String[] totalCoverage = null;
	String[] averageCoverage = null;
	int autosomalCount = 0;
	int xCount = 0;
	int elseCount = 0;
	public ChrXCheck(){
	}
	public void setNormalize(boolean n){
		this.normalize = n;
	}
	public void setSkipNonExistentSamples(boolean b){
		this.skipNonExistentSamples = b;
	}
	public void normalizeFile(File src, File dest, File list, File used) throws FileNotFoundException, IOException{
		this.setNormalize(true);
		this.load(src, list);
		this.check();
		PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(dest)));
		// header
		pw.print("Target\ttotal_coverage\taverage_coverage");
		for(String name: this.genders.keySet()){
			pw.print("\t" + name + "_total_cvg");
			pw.print("\t" + name + "_mean_cvg");
			pw.print("\t" + name + "_granular_Q1");
			pw.print("\t" + name + "_granular_median");
			pw.print("\t" + name + "_granular_Q3");
			pw.print("\t" + name + "_%_above_15");
		}
		pw.println("");
		for(int i = 0; i<this.autosomalCount; i++){
			pw.print(this.target[i] + "\t" + this.totalCoverage[i] + "\t" + this.averageCoverage[i]);
			for(String name: this.genders.keySet()){
				SampleDepth data = this.depth.get(name);
				pw.print("\t" + data.auto_total[i]);
				pw.print("\t" + data.auto_mean[i]);
				pw.print("\t" + data.auto_gq1[i]);
				pw.print("\t" + data.auto_med[i]);
				pw.print("\t" + data.auto_gq3[i]);
				pw.print("\t" + data.auto_a15[i]);
			}
			pw.println("");
		}
		for(int i = 0; i<this.xCount; i++){
			pw.print(
				this.target[i+this.autosomalCount] + "\t" +
				this.totalCoverage[i + this.autosomalCount] + "\t" +
				this.averageCoverage[i + this.autosomalCount]);
			for(String name: this.genders.keySet()){
				SampleDepth data = this.depth.get(name);
				pw.print("\t" + data.x_total[i]);
				pw.print("\t" + data.x_mean[i]);
				pw.print("\t" + data.x_gq1[i]);
				pw.print("\t" + data.x_med[i]);
				pw.print("\t" + data.x_gq3[i]);
				pw.print("\t" + data.x_a15[i]);
			}
			pw.println("");
		}
		pw.close();

		pw = new PrintWriter(new BufferedWriter(new FileWriter(used)));
		for(String name: this.usedGenders.keySet()){
			UsedInfo info = this.usedGenders.get(name);
			pw.println(info.name + "\t" + info.gender + "\t" + info.fvalue + "\t" + info.nfvalue);
		}
		pw.close();
	}
	//
	public static double variance(double[] buf, int fold){
		// mean
		double mean = 0.0;
		for(double d: buf){
			mean += (d*fold);
		}
		mean = mean/buf.length;

		// variance
		double v = 0.0;
		for(double d: buf){
			d = d*fold;
			v += ((d-mean)*(d-mean));
		}
		return 1.0/(buf.length-1) * v;
	}
	public static double variance(double[] buf){
		return variance(buf, 1);
	}
	public static double variance_norm(double[] buf, double sX, double sAuto){
		double stdX = Math.sqrt(sX);
		double stdAuto = Math.sqrt(sAuto);
		double stddevRatio = Math.sqrt(sAuto)/Math.sqrt(sX);
		double mean = mean(buf);
		double v = 0.0;
		for(double d: buf){
			double diff = (d - mean)*stddevRatio;
			v += (diff*diff);
		}
		return 1.0/(buf.length-1) * v;
	}
	public void check(){
		for(String name: genders.keySet()){
			SampleDepth sample = depth.get(name);
			String gender = genders.get(name);
			UsedInfo info = new UsedInfo();
			info.name = name;
			info.gender = gender;
			usedGenders.put(name, info);
			double sAuto = variance(sample.auto_mean);
            if(sAuto == 0){
            	continue;
            }
			double sX    = variance(sample.x_mean);
			double fixSX = (gender.equals("M") || gender.equals("male"))?
				(this.normalize? variance_norm(sample.x_mean, sX, sAuto):variance(sample.x_mean, 2)) : sX;
			double f = sX/sAuto;
			double fixF = fixSX/sAuto;
			System.out.println(name + ": " + gender + "\t" + f +  "\t" + fixF);
			info.fvalue = f;
			info.nfvalue = fixF;
			if(this.normalize && f < 0.7){ // automatically force normalization without considering gender
				if(gender.equals("female") || gender.equals("unknown")){
					info.gender = "male";
					// XHMMLogger.println("DISCRIPANCY: " + name + " " + gender + " F=" + f);
					System.err.println("DISCRIPANCY: " + name + " " + gender + " F=" + f);
				}
				System.err.println("Normalizing X depth: " + name);
				normalizeX(Math.sqrt(sAuto)/Math.sqrt(sX), sample);
				sX = variance(sample.x_mean);
				fixSX = variance_norm(sample.x_mean, sX, sAuto);
			    fixF = fixSX/sAuto;
				info.nfvalue = fixF;
			}
			if( f >= 0.7 && (gender.equals("male") || gender.equals("unknown"))){
				info.gender = "female";
				// XHMMLogger.println("DISCRIPANCY: " + name + " " + gender + " F=" + f);
				System.err.println("DISCRIPANCY: " + name + " " + gender + " F=" + f);
			}
			if( f > 1.25 ){
				info.gender = "super-female";
				System.err.println("WARNING: " + name + " " + "super-female" + " F=" + f);
			}
		}
	}
	public void normalizeX(double ratio, SampleDepth sampleDepth){
		double mean = mean(sampleDepth.x_mean);
		double autoMean = mean(sampleDepth.auto_mean);
		for(int i = 0; i<sampleDepth.x_mean.length; i++){
			double x = sampleDepth.x_mean[i];
			double diff = (x - mean)*ratio;
			sampleDepth.x_mean[i] = autoMean + diff;
			if(autoMean + diff < 0){ // avoid minus value
				sampleDepth.x_mean[i] = 0.0;
			}
		}
	}
	public static double mean(double[] buf){
		double mean = 0.0;
		for(double d: buf){
			mean += d;
		}
		return mean/buf.length;
	}
	public void load(File depthFile, File listFile) throws FileNotFoundException, IOException{
		BufferedReader br = new BufferedReader(new FileReader(listFile));
		String raw = null;
		while(null != (raw = br.readLine())){
			raw = raw.trim();
			if(raw.length() == 0){
				System.err.println("An empty line was found in " + listFile.getName());
			}
			String[] line = raw.split("\t");
			if(line.length < 2){
				System.err.println("An insufficient line was found in " + listFile.getName() + ": " + raw);
				System.exit(1);
			}
			// System.err.println(line[0]);
			genders.put(line[0], line[1]);
		}
		br.close();
		br = new BufferedReader(new FileReader(depthFile));
		String[] header = br.readLine().split("\t");
		
		// prepare buffer
		this.autosomalCount = 0;
		this.xCount = 0;
		this.elseCount = 0;
		while(null != (raw = br.readLine())){
			String[] line = raw.split("\t");
			// [0] Target
			// [1] total_coverage
			// [2] average_coverage
			if(line[0].matches("^[1-9].*")){
				autosomalCount++;
			}else if(line[0].startsWith("X")){
				xCount++;
			}else {
				elseCount++;
			}
		}
		System.err.println("Autosome: " + autosomalCount + " lines, ChrX: " + xCount);
		br.close();
		for(String sampleId: genders.keySet()){
			depth.put(sampleId, new SampleDepth(autosomalCount, xCount, elseCount));
		}
		// parse data
		target =  new String[autosomalCount + xCount + elseCount];
		totalCoverage = new String[autosomalCount + xCount + elseCount];
		averageCoverage = new String[autosomalCount + xCount + elseCount];
		br = new BufferedReader(new FileReader(depthFile));
		header = br.readLine().split("\t");
		int lineNo = 0;
		while(null != (raw = br.readLine())){
			String[] line = raw.split("\t");
			target[lineNo] = line[0];
			totalCoverage[lineNo] = line[1];
			averageCoverage[lineNo] = line[2];
			lineNo++;
			for(int col = 3; col < line.length; col+=6){
				SampleDepth d = depth.get(header[col].replaceAll("_total_cvg", ""));
				if(d == null && skipNonExistentSamples){
					continue;
				}else if(d == null){
					System.err.println(skipNonExistentSamples);
					System.err.println("[ERROR] The sample '" + header[col].replaceAll("_total_cvg", "") +"' does not exist in the list " + listFile.getName());
					System.exit(-1);
				}
				d.add(line[0], line[col], line[col+1], line[col+2], line[col+3], line[col+4], line[col+5]);
			}
		}
		br.close();
	}
	public static void main(String[] argv){
		try {
			String sampleIntervalSummary = null;
			String normalizedIntervalSummaryOut = null;
			String registeredGenderInfo = null;
			String deducedGenderOut = null;
			if(argv.length < 8){
				System.err.println("Usage: java ChrXCheck [-s sample_interval_summary(input)] [-n normalized_sample_interval_summary_out(output)] [-r registered_gender_info(input)] [-d deduced_gender_out(output)]\n");
				System.exit(0);
			}else {
				for(int i = 0; i<argv.length-1; i++){
					if(argv[i].equals("-s")){
						sampleIntervalSummary = argv[i+1];
					}else if(argv[i].equals("-n")){
						normalizedIntervalSummaryOut = argv[i+1];
					}else if(argv[i].equals("-r")){
						registeredGenderInfo = argv[i+1];
					}else if(argv[i].equals("-d")){
						deducedGenderOut = argv[i+1];
					}
				}
			}
			if(sampleIntervalSummary == null){
				System.err.println("Use -s option for sample interval summary file");
				System.exit(-1);
			}
			if(normalizedIntervalSummaryOut == null){
				System.err.println("Use -n option for normalized interval summary file");
				System.exit(-1);
			}
			if(registeredGenderInfo == null){
				System.err.println("Use -r option for registered gender list file");
				System.exit(-1);
			}
			if(deducedGenderOut == null){
				System.err.println("Use -d option for deduced gender list file");
				System.exit(-1);
			}

			// normalizeFile(File src, File dest, File list, File used)
			ChrXCheck c = new ChrXCheck();
			c.setSkipNonExistentSamples(true);
			c.normalizeFile(new File(sampleIntervalSummary), new File(normalizedIntervalSummaryOut), new File(registeredGenderInfo), new File(deducedGenderOut));
		}catch(FileNotFoundException fe){
			System.err.println("File not found:" + fe.getMessage());
		}catch(Exception e){
			e.printStackTrace();
		}
	}
}
