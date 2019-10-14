extern crate bio_utils;
extern crate last_tiling;
fn main() -> std::io::Result<()> {
    let args: Vec<_> = std::env::args().collect();
    let maf: Vec<_> = bio_utils::maf::Reader::from_file(&args[1])?
        .records()
        .filter_map(|e| e.ok())
        .collect();
    let tab: Vec<_> = last_tiling::parse_tab_file(&args[2])?;
    eprintln!("{}", maf.len());
    eprintln!("{}", tab.len());
    // Sanitization check.
    println!("Pm\tPSi\tPSd\tPEi\tPEd\tPid\tPdi");
    {
        let (p_m, p_s_i, p_s_d, p_e_i, p_e_d, p_i_d, p_d_i) = summarize_maf(tab);
        println!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            p_m, p_s_i, p_s_d, p_e_i, p_e_d, p_i_d, p_d_i
        );
    }
    {
        let (p_m, p_s_i, p_s_d, p_e_i, p_e_d, p_i_d, p_d_i) = summarize_tab(maf);
        println!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            p_m, p_s_i, p_s_d, p_e_i, p_e_d, p_i_d, p_d_i
        );
    }

    Ok(())
}
