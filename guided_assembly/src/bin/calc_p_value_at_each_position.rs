extern crate rust_htslib;
use rust_htslib::bam;
use rust_htslib::bam::Read;
fn main()->std::io::Result<()>{
    let args:Vec<_> = std::env::args().collect();
    // /3 is needed to convert 'multi-way substitution rate' into
    // one way error.
    let subst_rate:f64 = args[2].parse::<f64>().unwrap() / 3f64;
    let pileups = get_pileup(&args[1],subst_rate)?;
    println!("tid\tposition\tp_value\tdepth\tminor_allel_count");
    for (tid, position, p_value,depth,mac) in pileups{
        println!("{}\t{}\t{}\t{}\t{}",tid,position,p_value,depth,mac);
    }
    Ok(())
}

fn get_pileup(file:&str,subst_rate:f64)->std::io::Result<Vec<(u32,u32,f64,u32,u32)>>{
    let mut bam = bam::Reader::from_path(&std::path::Path::new(file))
        .map_err(|why|{
            eprintln!("{:?}",why);
            std::io::Error::new(std::io::ErrorKind::Other,"Error")
        })?;
    Ok(
        bam.pileup()
            .filter_map(|e|e.ok())
            .filter_map(|pileup| calc_p_value(&pileup,subst_rate))
            .collect()
    )
}
fn calc_p_value(pileup: &bam::pileup::Pileup,subst_rate:f64)->Option<(u32,u32,f64,u32,u32)>{
    let tid = pileup.tid();
    let pos = pileup.pos();
    let (depth, minor_count) = count(pileup)?;
    let mu = depth as f64 * subst_rate; // Expected number of minor allel count.
    if (minor_count as f64) < mu {
        Some((tid,pos,1.,depth,minor_count)) // If is not true, but keep it simple.
    }else{
        let delta_minus_one = minor_count as f64 / mu;
        let p_value = (mu.exp() / delta_minus_one.powf(delta_minus_one)).powf(mu);
        Some((tid, pos, p_value.min(1.),depth,minor_count))
    }
}
                                          
fn count(pileup:&bam::pileup::Pileup)->Option<(u32,u32)>{
    let _pos = pileup.pos();
    let _tid = pileup.tid();
    //eprintln!("Start {}\t{}",tid,pos);
    let (mut a, mut c, mut g, mut t) = (0,0,0,0);
    pileup.alignments()
        .filter(|align|!align.is_del() && !align.is_refskip())
        .filter(|align|{
            let flags = align.record().flags();
            flags & 0x800 == 0 && flags & 0x100 == 0
        })
        .filter_map(|align|{
            if align.record().seq().len()==0 || align.qpos().is_none(){
                // eprintln!("{}",align.record().flags());
                return None;
            }
            let base = align.record().seq()[align.qpos().unwrap()];
            // if align.record().is_reverse(){
            //     eprint!("({},-)", base as char);
            // }else{
            //     eprint!("({},*)",base as char);
            // }
            Some(base)
        })
        .for_each(|base| match base {
            base if base == b'a' || base == b'A' => a = a+1,
            base if base == b'c' || base == b'C' => c = c+1,
            base if base == b'g' || base == b'G' => g = g+1,
            base if base == b't' || base == b'T' => t = t+1,
            _ => {},
        });
    let depth = a + c + g + t;
    // eprintln!("\nA{},C{},G{},T{}",a,c,g,t);
    if depth == 0{
        None
    }else{
        Some((depth, second_largest(vec![a,c,g,t])))
    }
}

fn second_largest(data:Vec<u32>)->u32{
    let max = data.iter().max().unwrap();
    if let Some(max) = data.iter().filter(|&e| e < max).max(){
        *max
    }else{
        data[0]
    }
}
