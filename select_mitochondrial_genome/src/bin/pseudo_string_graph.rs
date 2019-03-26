use std::collections::HashMap;
#[derive(Clone,Debug)]
struct StringGraph<'a>{
    indexer:HashMap<String,usize>,
    edges:Vec<Vec<&'a str>>,
}

impl<'a> std::fmt::Display for StringGraph<'a>{
    fn fmt(&self,f:&mut std::fmt::Formatter) -> std::fmt::Result{
        writeln!(f, "There are {} nodes.",self.indexer.len())?;
        let numofedge:usize = self.edges.iter().map(|e|e.len()).sum();
        writeln!(f, "There are {} edges.",numofedge);
        writeln!(f, "ND Nodes. To See the nodes list, `grep ^ND`")?;
        for (key,val) in &self.indexer{
            writeln!(f, "ND {}\t{}",key,val)?;
        }
        writeln!(f, "ED Edges. Source -> Destination1, Destination2,.....")?;
        writeln!(f, "ED To See the edge list, `grep ^ED`")?;
        for (source,edges) in self.edges.iter().enumerate(){
            write!(f,"{} -> ",source)?;
            for dest in edges {
                let dest = self.indexer[*dest];
                write!(f,"{},",dest)?;
            }
            writeln!(f,"\x08")?;
        }
        Ok(())
    }
}

impl<'a> StringGraph<'a>{
    fn new(paf:String)->Self{
        StringGraph{
            indexer:HashMap::new(),
            edges:vec![],
        }
    }
}

fn test_paf()->String{
    let paf =
        "test\t100\t0\t50\t+\ttest2\t100\t50\t100\t25\t55\t60\n";
    paf.to_string()
}

#[test]
fn test(){
    let paf = test_paf();
    let string_graph = StringGraph::new(paf);
    println!("{}",string_graph);
}

fn open_paf(file:&str)->std::io::Result<String>{
    use std::io::Read;
    let mut input = String::new();
    let mut file = std::fs::File::open(&std::path::Path::new(file))?;
    file.read_to_string(&mut input)?;
    Ok(input)
}

fn main()->std::io::Result<()>{
    let args:Vec<_> = std::env::args().collect();
    let paf:String = open_paf(&args[1])?;
    let string_graph = StringGraph::new(paf);
    println!("{}",string_graph);
    println!("{}",test_paf());
    Ok(())
}
