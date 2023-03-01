use std::fs::File;
use std::io::Read;
use std::env;
use std::io::Write;
use bio::io::fasta;
use bincode::deserialize;
use std::time::Instant;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use std::collections::HashMap;

#[derive(Serialize, Deserialize)]
struct BuildSAResults {
    sequence:Vec<u8>,
    suffix_array: Vec<usize>,
    prefix_table: BTreeMap<Vec<u8>, Vec<u32>>,
}


fn main() -> std::io::Result<()> {
    let before = Instant::now();
    let args: Vec<String> = env::args().collect();
    let index: &str = &args[1];
    let query: &str = &args[2];
    let query_mode: Option<&str> = Some(&args[3]);
    let output= &args[4];
    let results = querysa(index,query,query_mode);
    println!("Elapsed time: {:.2?}", before.elapsed());
    //println!("{:?}",results);
    let mut file = File::create(output)?;
    for (seqid, hits) in results {
        //let hits_tab = hits.join("\t");
        let mut line = String::new();
        line.push_str(&seqid);
        line.push_str("\t");

        line.push_str(&hits.len().to_string());
        line.push_str("\t");

        for hit in hits {
            line.push_str(&hit.to_string());
            line.push_str("\t");
        }

        file.write_all(line.trim().as_bytes())?;
        file.write_all(b"\n")?;

        
    }

    Ok(())

}


fn querysa(index: &str,query:&str,query_mode: Option<&str>) -> HashMap<std::string::String, Vec<usize>> {
    let my_data = read_data(index);
    let sa_results = my_data.unwrap();
    let query_seqs = getquery(query);
    //let mut hits_vec = Vec::new();
    //let mut seqid_vec= Vec::new();
    let mut results = HashMap::new();
    

    match query_mode {
        Some("naive") => {
            //let mut hits_vec = Vec::new();
            //let mut seqid_vec= Vec::new();
            for (k,i) in query_seqs{
                //println!("Sequence name:{:?}",k);
                //seqid_vec.push(k.to_owned());
                match naive_binary_search(&sa_results,i.clone()) {
                    Some(hits) => {
                        //hits_vec.push(hits.clone());
                        results.insert(k.to_owned(), hits.clone());},
                    None => {results.insert(k.to_owned(), Vec::new());},//println!("Sequece:{:?},No hits",k),
                } 
                
            }
        },
        Some("simpaccel") => {
            //let mut hits_vec = Vec::new();
            //let mut seqid_vec= Vec::new();
            for (k,v) in query_seqs{
                //seqid_vec.push(k.to_owned());
                let hits =  acc_search(&sa_results,v.clone());
                //hits_vec.push(hits.clone());
                results.insert(k.to_owned(), hits.clone());
            }
        },
        None => println!("No query mode specified, performing default operation"),
        _ => println!("Invalid query mode specified"),
    }
    
    results
    
}






fn acc_search(sa_results:&BuildSAResults, query:Vec<u8>)-> Vec<usize>{
    
    let mut interval=Vec::new();
    match get_sa_interval(&sa_results,query.clone()) {
        Some(sa_interval) => {interval.push(sa_interval);},
        None => {},//println!("No matched suffix interval"),
    }
    
    let interval_tmp: Vec<u32> = interval
        .iter()                         
        .flat_map(|inner_vec| inner_vec.iter())  
        .cloned()                       
        .collect();     
    let suffix_array_interval: Vec<usize> = interval_tmp.iter().map(|&x| x as usize).collect();   

    let text=&sa_results.sequence;
    //let suffix_array = &sa_results.suffix_array; 
    //let mut hits = Vec::new();  
    let mut hits = Vec::new();  
    if suffix_array_interval.len() > 0{
        for i in suffix_array_interval{
            if text[i..].starts_with(query.as_slice()){
                hits.push(i);
            }

        }
        
    }
    hits     
    
}

fn read_data(index: &str) -> Result<BuildSAResults, Box<dyn std::error::Error>> {
    let mut file = File::open(index)?;
    let mut buf = Vec::new();
    file.read_to_end(&mut buf)?;
    let results: BuildSAResults = deserialize(&buf)?;
    Ok(results)
} 


fn getquery(query:&str) -> HashMap<std::string::String, Vec<u8>> {
    let reader = fasta::Reader::from_file(query).unwrap();
    let mut sequences = HashMap::new();
    for record in reader.records() {
        let record = record.expect("Error reading record");
        let sequence = record.seq();
        let id = record.id();
        sequences.insert(id.to_owned(),sequence.to_owned());
    }
    sequences
}


fn get_sa_interval(sa_results:&BuildSAResults,query:Vec<u8>) -> Option<&Vec<u32>> {
    //let suffix_array = &sa_results.suffix_array;
    let preftab = &sa_results.prefix_table;
    let keys: Vec<Vec<u8>> = preftab.keys().cloned().collect();
    let k = keys[0].len();
    let query_k = &query[..k];
    
    for key in keys {
        if query_k == &key {
            //return Some(key);
            let sa_interval = preftab.get(&key).unwrap();
            return Some(sa_interval);
        } 
    }  
    None
}





fn naive_binary_search(sa_results:&BuildSAResults, query:Vec<u8>)-> Option<Vec<usize>> {

    let suffix_array = &sa_results.suffix_array;
    let text=&sa_results.sequence;
    
    let mut left = 0;
    let mut right = suffix_array.len() - 1;
    let mut hit = Vec::new();

    while left <= right {
        let mid = (left + right) / 2;
        let suffix_start = suffix_array[mid];

        // Compare target string with the suffix starting at suffix_start.
        let suffix_end = suffix_start + query.len();
        
        let mut _suffix = Vec::new();
        if suffix_end < suffix_array.len() {
            _suffix = (&text[suffix_start..suffix_end]).to_vec();
        }else {
            _suffix = (&text[suffix_start..suffix_array.len()]).to_vec();
        }

        if query.as_slice() > &_suffix {
            left = mid + 1;
        } else if query.as_slice() < &_suffix {
            right = mid - 1;
        } else {
            hit.push(suffix_start);
            //return Some(suffix_start);
            let mut start = mid;
            let mut end = mid;
            //let mut hit = Vec::new();

            while start > 0 && text[suffix_array[start - 1]..].starts_with(query.as_slice()) {
                start -= 1;
                hit.push(suffix_array[start]);
            }

            while end < suffix_array.len() - 1 && text[suffix_array[end + 1]..].starts_with(query.as_slice()) {
                end += 1;
                hit.push(suffix_array[end]);
            }
            
            return Some(hit);
            
        }
        
    }
    
    None
}

