import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map.Entry;

public class Ex1 {
	private static HashMap<String, WordFrequency > tagWordFrequencies = new HashMap<String, WordFrequency >();
	private WordFrequency totals = new WordFrequency();
	
	public static void main(String[] args) throws IOException {
		new Ex1().run();
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
                    }
                }
            }
        }
	}
	
	private void run() throws IOException {
		init();
		
		
		System.out.println("Testing...");
		File g = new File("data/gene.test.p1");
		if(g.exists()) g.delete();
		g.createNewFile();
		BufferedWriter w = new BufferedWriter(new FileWriter(g));
		
		File f = new File("data/gene.test");
		BufferedReader r = new BufferedReader(new FileReader(f));
		String line;
		while((line = r.readLine()) != null) {
			if(! line.isEmpty()) {
				w.write(line + " " + getMLTag(line));
			}
			w.newLine();
		}
		r.close();
		w.close();
		
		System.out.println("Complete.");
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
