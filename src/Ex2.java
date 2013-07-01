import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.Stack;
import java.util.Map.Entry;

public class Ex2 {
	private static HashMap<String, WordFrequency > tagWordFrequencies = new HashMap<String, WordFrequency >();
	private static HashMap<String, Integer> ngramCount = new HashMap<String, Integer>();
	private static ArrayList<HashMap<String, Double>> viterbiProb = new ArrayList<HashMap<String, Double>>();
	private static ArrayList<HashMap<String, String>> viterbiTag = new ArrayList<HashMap<String, String>>();
	private static HashMap<String, Boolean> commonWords = new HashMap<String, Boolean>();
	
	private WordFrequency totals = new WordFrequency();
	
	public static void main(String[] args) throws IOException {
		new Ex2().run();
	}
	
	public void init() throws IOException {
        File f = new File("data/gene.counts");
        BufferedReader r = new BufferedReader(new FileReader(f));
        
        System.out.println("Training...");
        String line;
        while((line = r.readLine()) != null) {
            String[] parts = line.split(" ");
            if(parts[1].equals("WORDTAG")) {
                if(! tagWordFrequencies.containsKey(parts[2])) {
                    tagWordFrequencies.put(parts[2], new WordFrequency());
                }
                tagWordFrequencies.get(parts[2]).add(parts[3], Integer.parseInt(parts[0]));
                totals.add(parts[3], Integer.parseInt(parts[0]));
                commonWords.put(parts[3], true);
            } else if(parts[1].equals("1-GRAM")) {
                ngramCount.put(parts[2], Integer.parseInt(parts[0]));
            } else if(parts[1].equals("2-GRAM")) {
                ngramCount.put(parts[2] + " " + parts[3], Integer.parseInt(parts[0]));
            } else if(parts[1].equals("3-GRAM")) {
                ngramCount.put(parts[2] + " " + parts[3] + " " + parts[4], Integer.parseInt(parts[0]));
            }
        }
        r.close();
        
        System.out.println("Identifying Rare Words...");
        for(Entry<String, Integer> totalsEntry : totals.entrySet()) {
            if(totalsEntry.getValue() < 5) {
                for(Entry<String, WordFrequency> tagEntry : tagWordFrequencies.entrySet()) {
                    if(tagEntry.getValue().containsKey(totalsEntry.getKey())) {
                        int t = tagEntry.getValue().get(totalsEntry.getKey());
                        tagEntry.getValue().remove(totalsEntry.getKey());
                        tagEntry.getValue().add("_RARE_", t);
                        commonWords.remove(totalsEntry.getKey());
                    }
                }
            }
        }
	}
	
	private void run() throws IOException {
		init();
		
		System.out.println("Testing...");
		File g = new File("data/gene.test.p2");
		if(g.exists()) g.delete();
		g.createNewFile();
		BufferedWriter w = new BufferedWriter(new FileWriter(g));
		
		File f = new File("data/gene.test");
		BufferedReader r = new BufferedReader(new FileReader(f));
		String line;
		Stack<String> sentenceWords = new Stack<String>();
		
		int i = 0;
		while((line = r.readLine()) != null) {
		    sentenceWords.push(line);
		    runViterbi(++i, line);
		    
		    if(line.isEmpty()) {
		        sentenceWords.pop();
		        System.out.println("\nEnd of sentence. Finding best tag sequence:");
		        String sentenceOutput = "";
		        String endTag = "STOP";
		        while(i > 1) {
		            //System.out.println(sentenceWords.peek());
		            //System.out.println("i = " + i);
		            double pBest = 0;
		            String tBest = null;
    		        for(Entry<String, Double> entry : viterbiProb.get(i-1).entrySet()) {
    		            String[] parts = entry.getKey().split(" ");
    		            if(parts[1].equals(endTag)) {
    	                    //System.out.println("\t" + parts[0] + " " + parts[1] + ", p = " + entry.getValue());
    		                if(entry.getValue() > pBest) {
    		                    pBest = entry.getValue();
    		                    tBest = parts[0];  
    		                }  
    		            }
    		        }
    		        endTag = tBest;
    		        sentenceOutput = sentenceWords.pop() + " " + tBest + "\n" + sentenceOutput;
    		        //System.out.println("... " + sentenceTags + "\n");
    		        i--;
		        }
		        i=0;
		        System.out.println("\n" + sentenceOutput);
		        w.write(sentenceOutput);
		        w.newLine();
		        //break;
		    }
		    
			/*if(! line.isEmpty()) {
				w.write(line + " " + getMLTag(line));
			}
			w.newLine();*/
		}
		r.close();
		w.close();
		
		System.out.println("Complete.");
	}
	
	private void runViterbi(int pos, String word) {
	   // System.out.println("Viterbi algorithm, position " + pos + ". Word = " + word);
	    HashMap<String, Double> probMap = new HashMap<String, Double>();
	    viterbiProb.add(pos-1, probMap);
	   // HashMap<String, String> tagMap = new HashMap<String, String>();
       // viterbiTag.add(pos-1, tagMap);
	    
        HashMap<String, Double> prevProbMap = null;
        if(pos > 1) {
            prevProbMap = viterbiProb.get(pos-2);
        }
        
	    String[] tags = {"O", "I-GENE"};
	    String[] startTags = {"*"};
	    String[] stopTags = {"STOP"};
	    String[] tags1 = (pos > 2) ? tags : startTags;
	    String[] tags2 = (pos > 1) ? tags : startTags;
	    String[] tags3 = word.isEmpty() ? stopTags : tags;
	    
	    for(String tag3 : tags3) {
	        for(String tag2 : tags2) {
	            double pBest = 0;
	            String tBest = null;
	            for(String tag1 : tags1) {
	                double q = getTagProbabilityFromTrigram(tag1, tag2, tag3);
	                double r = (prevProbMap == null) ? 1 : prevProbMap.get(tag1 + " " + tag2);
	                double e;
	                if(word.isEmpty()) {
	                    e = 1;
	                } else {
    	                WordFrequency twf = tagWordFrequencies.get(tag3);
    	                if(! commonWords.containsKey(word)) {
    	                    //System.out.println(word + " is RARE!");
    	                    //System.out.println(tag3 + " has rare count:" + twf.get("_RARE_"));
    	                    e = (double)twf.get("_RARE_") / twf.total;
    	                } else if(twf.containsKey(word)) {
    	                    e = (double)twf.get(word) / twf.total;
    	                } else {
    	                    e = 0;
    	                }
	                }
	                double p = r * q * e;
	                if(pBest < p) {
	                    pBest = p;
	                    tBest = tag1;
	                }
	                //System.out.println("\t" + tag1 + " " + tag2 + " [" + pos + "] " + tag3 + "(" + word + ") has p = " + p + " ( = " + r + " * " + q + " * " + e + ")");
	            }
	            probMap.put(tag2 + " " + tag3, pBest);
	            //tagMap.put(tag2 + " " + tag3, tBest);
	            //System.out.println("\t ___ " + tag2 + " [" + pos + "] " + tag3 + "(" + word + ") has best p = " + pBest + " with tag " + tBest + "\n");
	        }
	    }
	}
	
	
	private double getTagProbabilityFromTrigram(String tag1, String tag2, String tag) {
	    if(ngramCount.containsKey(tag1 + " " + tag2 + " " + tag)) {
	        double c3 = ngramCount.get(tag1 + " " + tag2 + " " + tag); 
	        double c2 = ngramCount.get(tag1 + " " + tag2);
	        return (c3/c2);
	    }
	    return 0;
	}
	
	private String getMLTag(String word) {
		double eBest = 0;
		String tBest = null;
		for(Entry<String, WordFrequency> entry : tagWordFrequencies.entrySet()) {
			if(entry.getValue().containsKey(word)) {
				double e = (double)entry.getValue().get(word) / entry.getValue().getTotal();
				if(e > eBest) {
					eBest = e;
					tBest = entry.getKey();
				}
			}
		}
		
		if(tBest == null) {
			return getMLTag("_RARE_");
		}
		return tBest;
	}
	
	
	private class WordFrequency extends HashMap<String, Integer> {
		private int total = 0;
		
		public int getTotal() {
			return total;
		}
		
		@Override
		public Integer put(String key, Integer value) {
			total += value;
			return super.put(key, value);
		}
			
		public Integer add(String key, Integer value) {
			total += value;
			if(! containsKey(key)) {
				return super.put(key, value);
			} else {
				return super.put(key, value + get(key));
			}
		}
		
		@Override
		public Integer remove(Object key) {
			if(containsKey(key)) {
				total -= get(key);
			}
			return super.remove(key);
		}
		
	}

}
