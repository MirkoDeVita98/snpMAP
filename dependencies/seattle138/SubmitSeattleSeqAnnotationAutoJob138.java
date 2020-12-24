// SubmitSeattleSeqAnnotationAutoJob138.java
// Peggy Robertson, Nickerson Lab, Dept. of Genome Sciences, University of Washington
// 8/2016
// example client program for querying SeattleSeq Annotation from a local computer; runs with Java 1.5 or higher

/*
 
steps for using:
  
1. Modify this code to put your own inputFile name, outputFile name, eMail, compression choice, and timeout minutes in.
   The code will not compile until an email address is added.
   In the submitTheInputFile function, choose annotation parameters.

2. Execute "java -version" to make sure java is installed and has a version of at least 1.5 (plus any subversion).
 
3. Acquire some jar files if not present (newer versions may be ok):
	httpunit.jar
	nekohtml-0.9.5 or jtidy-4aug2000r7-dev.jar
	js-1.6R5.jar
	xercesImpl-2.6.1.jar

4. Set your classpath variable with your own path substituted for "path" below (this is for bash, and can be put in your .bashrc file so you only have to do it once):
export CLASSPATH=./:/Users/sat/Desktop/ETH/ETH_PROJECTS/comp_bio/team02/P2/seattle138/httpunit-1.7.jar:/Users/sat/Desktop/ETH/ETH_PROJECTS/comp_bio/team02/P2/seattle138/js-1.6R5.jar:/Users/sat/Desktop/ETH/ETH_PROJECTS/comp_bio/team02/P2/seattle138/nekohtml-0.9.5.jar:/Users/sat/Desktop/ETH/ETH_PROJECTS/comp_bio/team02/P2/seattle138/xercesImpl-2.6.1.jar
export CLASSPATH=./:/home/pierobartolo/ETH/comp_bio/team02/P2/dependencies/seattle138/httpunit-1.7.jar:/home/pierobartolo/ETH/comp_bio/team02/P2/dependencies/seattle138/js-1.6R5.jar:/home/pierobartolo/ETH/comp_bio/team02/P2/dependencies/seattle138/nekohtml-0.9.5.jar:/home/pierobartolo/ETH/comp_bio/team02/P2/dependencies/seattle138/xercesImpl-2.6.1.jar

5. Compile:
javac SubmitSeattleSeqAnnotationAutoJob138.java
 
6. Execute:
java SubmitSeattleSeqAnnotationAutoJob138
 
If your job fails and you see a line like "downloadURL = <html>", you may not have chosen the right input format.
 
*/

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.BufferedWriter;
import java.io.PrintWriter;
import java.io.FileWriter;
import java.io.FileOutputStream;

import java.util.ArrayList;
import java.util.zip.GZIPInputStream;

import com.meterware.httpunit.*;
import com.meterware.httpunit.protocol.*;

public class SubmitSeattleSeqAnnotationAutoJob138 {

	
	// ----------------------------- the parameters in this section must be customized

	int timeoutMinutes = 120;	// increment this for long jobs
	private String eMail ="pdebartol@ethz.ch" ;	// must be filled in; an email will not normally be sent, but error messages will be
	
	// examples for no output compression
	public String inputFilename = "/home/pierobartolo/ETH/comp_bio/team02/P2/cache/small_sSNP.csv";	// local file to be submitted
 	public String outputFilename = "/home/pierobartolo/ETH/comp_bio/team02/P2/cache/test_sSNP.seattle38.clinvar.vcf";	// local file where results will be written)
	private String institutionalString = "UW";	// an identifier that will be used by our server for writing files (use only alpanumerical characters)
	private boolean useCompressionForOutput = false;

	// examples for output compression
	//private String inputFilename = "/Users/pdr/solexa/dataForFunctionAnalysis/stdout9.alternateFormat.hg19.vcf";	// local file to be submitted
 	//private String outputFilename = "/Users/pdr/misc/stdout9.alternateFormat.hg19.results.compress.txt.gz";	// local file where results will be written)
	//private String institutionalString = "UW";	// an identifier that will be used by our server for writing files (use only alpanumerical characters)
	//private boolean useCompressionForOutput = true;

	// ----------------------------- 

	
	final int BUFFER_SIZE = 65536;

	private String URL = "http://snp.gs.washington.edu/SeattleSeqAnnotation138/";	// SeattleSeq Annotation URL

	private String downloadURL = "";	// URL to get result file
	private String progressURL = "";	// URL to monitor when job is done

	WebConversation wc = null;

	public SubmitSeattleSeqAnnotationAutoJob138() {
	}
	
	public void runJob() {
		try {
			// create the conversation object
			wc = new WebConversation();

 			// obtain the main page on the SeattleSeq Annotation site
			System.out.println("ask for the main page");
			WebResponse response = getMainPageResponse(wc);
			
			// submit the request
			System.out.println("submit the file");
			response = submitTheInputFile(response, inputFilename, wc);	// fill out the form on the main page with the file, the email address, and the options
			
			// set the URLs for download and progress
			boolean isValid = setURLsFromSource(response);
			if (!isValid) {
				System.out.println("The attempt to get URLs for monitoring progress failed.  Most likely your submitted file does not have an autoFile line in it, or is not of the expected format.");
			}
			else if (downloadURL.contains("ABORT") || progressURL.contains("ABORT")) {
				System.out.println("downloadURL = " + downloadURL + "\nprogressURL = " + progressURL);
				System.out.println("ERROR: The SeattleSeqAnnotation file was aborted.  The job may be too large or too many simultaneous jobs may have been submitted.");
			}
			else {	// the normal case
				// e.g.
				// http://gvsbatch.gs.washington.edu/SeattleSeqAnnotation138/BatchFileDownloadServlet?file=testAuto.123456789.txt&download=plainText
				// http://gvsbatch.gs.washington.edu/SeattleSeqAnnotation138/BatchProgressServlet?progressFile=progress.testAuto.123456789.txt&auto=true
				System.out.println("downloadURL = " + downloadURL + "\nprogressURL = " + progressURL);
				
				// wait for the job to finish, monitoring with the progress file
				boolean isSuccessfulFinish = waitForJobToFinish();	// use progressURL to detect job completion
				
				// download the result file
				if (isSuccessfulFinish) {
					WebRequest requestDownload = new GetMethodWebRequest(downloadURL);	// SeattleSeq Annotation uses http GET here
					WebResponse responseDownload = wc.getResponse(requestDownload);
					writeSourceToFile(responseDownload);	// write the downloaded results to a local file
					System.out.println("download complete");
				}
				else {
					System.out.println("PROBLEM: could not detect successful finish of job");
				}
			}
		} catch (Exception ex) {
			System.err.println("ERROR runJob - " + ex);
			ex.printStackTrace();
		}
	}

	private WebResponse getMainPageResponse(WebConversation wc) {
		// obtain the main page on the SeattleSeq Annotation site
		WebResponse responseResult = null;
		try {
			WebRequest request = new GetMethodWebRequest(URL);
			responseResult = wc.getResponse(request);
		} catch (Exception ex) {
			System.err.println("ERROR getMainPageResponse - " + ex);
		}	
		return responseResult;
	}		

	private WebResponse submitTheInputFile(WebResponse response, String inputFilename, WebConversation wc) {
		WebResponse responseInput = null;
		try {
			// fill out the form for submitting a file
			WebForm form = response.getFormWithName("GenotypeSummary");				// GenotypeSummary is the name of the form on SeattleSeq Annotation
			SubmitButton button = form.getSubmitButton("gFetch");					// gFetch is the name of the submit button			
			UploadFileSpec fileSpec = new UploadFileSpec(new File(inputFilename));
			
			WebRequest webRequest = form.newUnvalidatedRequest(button);
			webRequest.setParameter("GenotypeFile", new UploadFileSpec[] { fileSpec });	// GenotypeFile is the name of the file field
			webRequest.setParameter("EMail", eMail);										// EMail is the name of the email field
			// set various options; more parameter names may be found by viewing the html source for the home page
			//webRequest.setParameter("fileFormat", "Maq");	// use this to choose input format Maq
			webRequest.setParameter("fileFormat", "VCFSNVsAndIndels");                        
			webRequest.setParameter("autoFile", institutionalString);
			// webRequest.setParameter("HapMapFreqType", "HapMapFreqRef");	// use this to choose HapMap frequencies by reference-allele
            
            if (useCompressionForOutput) {
                webRequest.setParameter("compressAuto", "true");
            }

			// for debugging
			//String[] parameterNames = webRequest.getRequestParameterNames();
			//for (String name : parameterNames) {
			//	System.out.println(name + " " + webRequest.getParameterValues(name));
			//}
			
			ArrayList <String> columnsArray = new ArrayList <String> ();	// add any additional annotation columns wanted
			columnsArray.add("sampleAlleles");
			//columnsArray.add("allelesDBSNP");
			//columnsArray.add("scorePhastCons");
			//columnsArray.add("consScoreGERP");
            
            // CADD scores are Copyright 2013 University of Washington and Hudson-Alpha Institute for Biotechnology (all rights reserved)
            // but are freely available for all academic, non-commercial applications. For commercial licensing information contact
            // Jennifer McCullar (mccullaj@uw.edu). CADD is currently developed by Martin Kircher, Daniela M. Witten, Gregory M. Cooper,
            // and Jay Shendure ( http://cadd.gs.washington.edu/ ).
			//columnsArray.add("scoreCADD");
            
			//columnsArray.add("CNV");
			//columnsArray.add("dbSNPValidation");
			columnsArray.add("distanceToSplice");
			//columnsArray.add("microRNAs");
			columnsArray.add("tfbs");
			String [] columns = columnsArray.toArray(new String[columnsArray.size()]);	
			
			webRequest.setParameter("columns", columns);
			
			// click the submit button
			responseInput = wc.getResponse(webRequest);
		} catch (Exception ex) {
			System.err.println("ERROR submitTheInputFile - " + ex);
		}
		
		return responseInput;
	}

	private void writeSourceToFile(WebResponse response) {
		if (useCompressionForOutput) {
			System.out.println("The returned file is expected to be gzipped.");
			writeCompressedSourceToFile(response);
		}
		else {
			System.out.println("The returned file is expected to be plain text.");
			writeUncompressedSourceToFile(response);
		}
	}
	
	private void writeUncompressedSourceToFile(WebResponse response) {
		// write the results to a local file
		try {
			PrintWriter writer = new PrintWriter(new BufferedWriter(new FileWriter(outputFilename), BUFFER_SIZE));	// writes to local file
			InputStream inputStream = response.getInputStream();
			InputStreamReader reader = new InputStreamReader(inputStream);
			BufferedReader bufferedReader = new BufferedReader(reader);	// reads the plain-text source from the server
			String line1 = null;
			while((line1 = bufferedReader.readLine()) != null) {	// read a line, then write a line
				writer.println(line1);
			}
			inputStream.close();
			bufferedReader.close();
			writer.flush();
			writer.close();
		} catch (Exception ex) {
			System.err.println("ERROR writeUncompressedSourceToFile - " + ex);
			ex.printStackTrace();
		}		
	}

	private void writeCompressedSourceToFile(WebResponse response) {
		// write the results to a local file
		try {
			String contentType = response.getContentType();
			InputStream inputStream = response.getInputStream();
			if (contentType.equals("application/x-gzip") || contentType.equals("application/x-download")) {
				String compressedOutputFilename = outputFilename;
				if (!compressedOutputFilename.endsWith(".gz")) {
					compressedOutputFilename += ".gz";
				}
				File outputFile = new File(compressedOutputFilename);
				FileOutputStream outputStream = new FileOutputStream(outputFile);
				byte[  ] buf = new byte[4 * 1024];  // 4K char buffer
				int bytesRead;
				while ((bytesRead = inputStream.read(buf)) != -1) {
					outputStream.write(buf, 0, bytesRead);
				}
				outputStream.flush();
				outputStream.close();
			}
			else {
				System.err.println("ERROR writeCompressedSourceToFile: a compressed content type for the returned file was not detected; check the # compress line of the input file");
			}
		} catch (Exception ex) {
			System.err.println("ERROR writeCompressedSourceToFile - " + ex);
			ex.printStackTrace();
		}		
	}
	
	private boolean setURLsFromSource(WebResponse response) {
		// process the text from the submission-acknowledge http message, and extract the links for monitoring progress and downloading the result file
		boolean isValid = true;
		try {
			InputStream inputStream = response.getInputStream();
			InputStreamReader reader = new InputStreamReader(inputStream);
			BufferedReader bufferedReader = new BufferedReader(reader);
			String line1 = bufferedReader.readLine();
			String[] parts = line1.split(",");	// the first line in the source contains the download and progress URLs, comma-separated
			if (parts.length < 2 && line1.contains("html")) {
				isValid = false;
			}
			else {	// the normal case
				downloadURL = parts[0];
				progressURL = parts[1];
			}
			inputStream.close();
			bufferedReader.close();
		} catch (Exception ex) {
			System.err.println("ERROR setURLsFromSource - " + ex);
			ex.printStackTrace();
		}
		return isValid;
	}
	
	private boolean waitForJobToFinish() {
		// keep reading the progress file until the job is done, or it times out
		// the web response is two comma-separated numbers; the first is the number of lines completed, the second is the total number of lines
		// at the beginning, before any lines are processed, the response is 0,0

		boolean isSuccessfulFinish = false;
		int sleepMilliseconds = 10000;	// 10 seconds; this could be increased for long jobs
		int maxCycles = timeoutMinutes * 60000 / sleepMilliseconds;
		
		try {
			System.out.println("begin wait");
			int numberTries = 0;
			
			while (true) {
				numberTries++;
				if (numberTries > maxCycles) {	// failed to finish in alloted time
					System.out.println("ERROR: the annotation timed out after " + numberTries + " cycles, time-out minutes is " + timeoutMinutes);
					break;
				}
				
				// read information from the progress file on the server
				WebRequest requestProgress = new GetMethodWebRequest(progressURL);	//BatchProgressServlet uses GET
				WebResponse responseProgress = wc.getResponse(requestProgress);
				InputStreamReader reader = new InputStreamReader(responseProgress.getInputStream());
				BufferedReader bufferedReader = new BufferedReader(reader);
				String line1 = bufferedReader.readLine();
				String[] parts = line1.split(",");
				int numberLinesFinished = Integer.valueOf(parts[0]).intValue();
				int numberLinesTotal = Integer.valueOf(parts[1]).intValue();
				bufferedReader.close();
				System.out.println("cycle " + numberTries + " lines finished " + numberLinesFinished + " out of total " + numberLinesTotal);
				
				// see if finished	
				if (numberLinesFinished == numberLinesTotal && numberLinesFinished != 0) {	
					isSuccessfulFinish = true;
					break;
				}
				Thread.sleep(sleepMilliseconds);
			}
			Thread.sleep(1500);	// give it 5 more seconds, to be sure
			System.out.println("end wait");
		} catch (Exception ex) {
			System.err.println("ERROR waitForJobToFinish - " + ex);
		}
		return isSuccessfulFinish;	
	}
	
	public static void main(String args[]) throws Exception {
		SubmitSeattleSeqAnnotationAutoJob138 myself = new SubmitSeattleSeqAnnotationAutoJob138();
		myself.inputFilename = args[0];
		myself.outputFilename = args[1];
		myself.runJob();
	}
}
