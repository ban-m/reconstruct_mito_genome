const width = 1000;
const height = 1000;
// Margin in radian
const theta_margin = 0.01;
// const theta_margin = 0;
const gap_position = 0.05;
// the height of a contigs.
const contig_thick = 10; 
const coverage_thick = 5;
const gap_jitters = d3.randomNormal(0,0.01);
const read_thick = 4;
const eplison = 5.001;
const jitters = d3.randomNormal(0,eplison);
const confluent_margin = 0.01;
// Radius
const contig_radius = 350;
const coverage_min = contig_radius+contig_thick;
const coverage_max = 450;
const handle_points_radius = 100;
const read_radius = contig_radius-30;
const gap_min_radius = read_radius;
const gap_max_radius = contig_radius-3;
const gap_min = 500;
const gap_max = 2000;
const gap_scale = d3.scaleLog()
      .domain([gap_min,gap_max])
      .range([gap_min_radius, gap_max_radius])
      .clamp(true);

// Circle radius
const min_radius = 1;
const max_radius = 8;
const min_read_num = 2;
const offset = 50;
const svg = d3.select("#plot")
      .append("svg")
      .attr("width",width)
      .attr("height",height);
const contigs_layer = svg.append("g")
      .attr("transform", `translate(${width/2},${height/2+offset})`)
      .attr("class","contigs");
const coverage_layer = svg.append("g")
      .attr("transform", `translate(${width/2},${height/2+offset})`)
      .attr("class","coverages");
const temp_coverage_layer = svg.append("g")
      .attr("transform", `translate(${width/2},${height/2+offset})`)
      .attr("class","temp-coverages");
const start_stop_layer = svg.append("g")
      .attr("transform", `translate(${width/2},${height/2+offset})`)
      .attr("class","start-stop-read");
const read_layer = svg.append("g")
      .attr("transform", `translate(${width/2},${height/2+offset})`)
      .attr("class","read");
const cr_layer = svg.append("g")
      .attr("transform", `translate(${width/2},${height/2+offset})`)
      .attr("class","critical-region");
const tooltip = d3.select("body").append("div")
      .attr("class", "tooltip")
      .style("opacity", 0);
const info = d3.select("#info");

const calcScale = (contigs) => {
    // Input: Array of JSON object
    // Output: d3 Scale object
    // Requirements: each input object should have "length" attribute
    // Convert base pair into radian
    const num = contigs.length;
    const total = contigs.map(c => c.length).reduce((x,y)=>x+y);
    return d3.scaleLinear()
        .domain([0,total])
        .range([0,2 * Math.PI - num * theta_margin]);
};

const calcStartPosition = (contigs)=>{
    // Input: Array of JSON object.
    // Output: Array[Num]
    // Requirements: each input object should have "length" attribute
    // Map from the index(id) of contig into the start position of the contig(in radian).
    const scale = calcScale(contigs);
    const max = Math.max(...contigs.map(c => c.id));
    let cum_sum = 0;
    let start_pos = new Array(max);
    for (const contig of contigs){
        start_pos[contig.id] = scale(cum_sum) + contig.id * theta_margin;
        cum_sum += contig.length;
    }
    return start_pos;
};

const calcHandlePoints = (start_pos) => {
    // Input: Array[Num]
    // Output: Array[Array[Num]]
    // Requirement: None
    // Map from combinations of ID to the handle points of them.
    let handle_points = new Array();
    for (let i = 0 ; i < start_pos.length; i ++){
        handle_points.push(new Array(start_pos.length));
    }
    const max = start_pos.length-1;
    start_pos.forEach((v1,k1)=>{
        start_pos.forEach((v2,k2)=>{
            const next1 = (k1 == max ? Math.PI * 2 - theta_margin : start_pos[k1]);
            const next2 = (k2 == max ? Math.PI * 2 - theta_margin : start_pos[k2]);
            const val = (next1 + next2 + v1 + v2)/4 - theta_margin/2;
            handle_points[k1][k2] = val;
        });
    });
    return handle_points;
};

const calcCovScale = (contigs)=>{
    // Input: Array on JSON object
    // Output: d3.scale object
    // Requirements: each input object should have "length" attribute
    // Scale for convert coverage into radius.
    const max = Math.max(...contigs.map(contig => Math.max(...contig.coverages)));
    // const min = Math.min(...contigs.map(contig => Math.min(...contig.coverages)));
    return d3.scaleLinear()
        .domain([0,max])
        .range([coverage_min,coverage_max]);
};

const calcReadNumScale = (contigs) => {
    // Input: Array on JSON object
    // Output: d3.scale object
    // Requirements: Each object in the argument should have an array of integer, which is
    // named "start_stop."
    // Calculate the scale for start/stop vizualization.
    const total = contigs.flatMap(c => c.start_stop).reduce((x,y) => x+y);
    const num = contigs.map(c => c.start_stop.length).reduce((x,y)=> x+y);
    const max = Math.max(...contigs.flatMap(c => c.start_stop));
    console.log("mean", total/ num);
    console.log("max", max);
    return d3.scaleLog()
        .domain([min_read_num,max])
        .range([min_radius,max_radius])
        .clamp(true);
};

const readToPath = (read,handle_points,bp_scale,start_pos,unit_length)=>{
    // Input: JSON object, Array[Array[Num]], d3.scale, Array[Num], Num
    // Output: String
    // Requirements: read should have units attribute, each of which elements
    // should have either "G"(for Gap) or "E"(for Encode)
    // Note that E consists of 
    let path = d3.path();
    let units = Array.from(read.units).reverse();
    const r = read_radius;
    let gap = 0;
    let unit = {};
    while(!unit.hasOwnProperty("E")){
        unit = units.pop();
        if (unit == undefined){
            return "";
        }else if (unit.hasOwnProperty("G")){
            gap = unit.G;
        }
    };
    // Current ID of the contig 
    let current_contig = unit.E[0];
    let start = start_pos[current_contig] - Math.PI/2;
    let s_radian = start + bp_scale(unit_length*unit.E[1]);
    let e_radian = start + bp_scale(unit_length*unit.E[3]);
    if (gap != 0){
        path.moveTo(gap_scale(gap) * Math.cos(s_radian), gap_scale(gap)*Math.sin(s_radian));
        path.lineTo(read_radius * Math.cos(s_radian), read_radius * Math.sin(s_radian));
    }else{
        path.moveTo(read_radius * Math.cos(s_radian), read_radius * Math.sin(s_radian));
    }
    path.arc(0,0,r,s_radian, e_radian);
    let current_unit = unit.E[3];
    for (unit of units.reverse()){
        if (unit.hasOwnProperty("G")){
            continue;
        }
        const diff = Math.abs(unit.E[1]-current_unit);
        current_unit = unit.E[1];
        if (unit.E[0] == current_contig && diff < 100){
            s_radian = start + bp_scale(unit_length*unit.E[1]);
            e_radian = start + bp_scale(unit_length*unit.E[3]);
            path.arc(0,0,r,s_radian,e_radian);
            current_unit = unit.E[3];
        }else{
            // Change contig. Connect them.
            start = start_pos[unit.E[0]] - Math.PI/2;
            s_radian = start + bp_scale(unit_length*unit.E[1]);
            e_radian = start + bp_scale(unit_length*unit.E[3]);
            // Bezier Curve to new point from here.
            const control_radius = handle_points[current_contig][unit.E[0]] - Math.PI/2;
            const control_x = handle_points_radius*Math.cos(control_radius);
            const control_y = handle_points_radius*Math.sin(control_radius);
            current_contig = unit.E[2];
            path.quadraticCurveTo(control_x,control_y,r*Math.cos(s_radian),r*Math.sin(s_radian));
            path.arc(0,0,r,s_radian,e_radian);
            current_unit = unit.E[3];
        }
    }
    return path.toString();
};

const selectRead = (read,unitlen) => {
    // Input: JSON object, Num
    // Output: boolean
    // Requirements: input object should have "units" property,
    // which is actually vector of object with "Gap" or "Encode" property.
    // Filter read as you like.
    // const from = 0;
    // const to = 1;
    // const set = new Set(read.units.filter(u => u.hasOwnProperty("E")).map(u => u.E[0]));
    // const max_gap = Math.max(...read.units.filter(u => u.hasOwnProperty("G")).map(u => u.G));
    // return read.cluster.includes(24);
    // return read.cluster.length == 0;
    return true;
};

const getNumOfGapRead = reads => {
    // Input: [JSON object]
    // Output: Num
    // Requirements: each element should be 'read' object.
    // Return numbers of reads which is just Gap.
    return reads.filter(read => {
        let units = Array.from(read.units);
        let unit = {};
        while(!unit.hasOwnProperty("E")){
            unit = units.pop();
            if (unit == undefined){
                return true;
            }
        };
        return false;
    }).length;
};



// Below, critical object is a json ob
// {'CP': {'contig1': {'contig': 0,
//    'start_unit': 132,
//    'end_unit': 500,
//    'direction': 'UpStream'},
//   'contig2': {'contig': 0,
//    'start_unit': 1223,
//    'end_unit': 2432,
//    'direction': 'DownStream'}}}
// {'CR': {'pos': {'contig': 0,
//    'start_unit': 132,
//    'end_unit': 500,
//    'direction': 'UpStream'}}}

const criticalpairToPath = (cp, handle_points, bp_scale,start_pos, unit_length)=>{
    const r = read_radius;
    let path = d3.path();
    // Move to contig1
    const contig1 = cp["contig1"];
    const contig1_start_angle = start_pos[contig1["contig"]] - Math.PI/2;
    const start_angle_1 = contig1_start_angle + bp_scale(unit_length*contig1["start_unit"]);
    const end_angle_1 = contig1_start_angle + bp_scale(unit_length*contig1["end_unit"]);
    path.moveTo(r * Math.cos(start_angle_1), r * Math.sin(start_angle_1));
    path.arc(0,0,r,start_angle_1, end_angle_1);
    // Bezier Curve to contig2.
    const contig2 = cp["contig2"];
    const contig2_start_angle = start_pos[contig2["contig"]] - Math.PI/2;
    const start_angle_2 = contig2_start_angle + bp_scale(unit_length*contig2["start_unit"]);
    const control_radius = handle_points[contig1["contig"]][contig2["contig"]] - Math.PI/2;
    const control_x = handle_points_radius*Math.cos(control_radius);
    const control_y = handle_points_radius*Math.sin(control_radius);
    path.quadraticCurveTo(control_x,control_y,r*Math.cos(start_angle_2),r*Math.sin(start_angle_2));
    const end_angle_2 = contig2_start_angle + bp_scale(unit_length*contig2["end_unit"]);
    path.arc(0,0,r, start_angle_2, end_angle_2);
    path.quadraticCurveTo(control_x,control_y,r*Math.cos(start_angle_1),r*Math.sin(start_angle_1));
    return path.toString();
};

const confluentregionToPath = (cr, handle_points, bp_scale,start_pos, unit_length)=>{
    const r = read_radius + 50;
    let path = d3.path();
    const contig = cr["pos"];
    const contig_start_angle = start_pos[contig["contig"]] - Math.PI/2;
    const start_angle = contig_start_angle + bp_scale(unit_length*contig["start_unit"]) - confluent_margin;
    const end_angle = contig_start_angle + bp_scale(unit_length*contig["end_unit"]) + confluent_margin;
    path.moveTo(r * Math.cos(start_angle), r * Math.sin(start_angle));
    path.lineTo(r * Math.cos(end_angle), r * Math.sin(end_angle));
    return path.toString();
};

const crToPath = (cr, handle_points, bp_scale,start_pos, unit_length)=>{
    // Input: JSON object, JSON object, Integer
    // Output: String
    // Requirements: Critical region object, scales
    // Return the path btw critical region, or confluent path.
    if (cr.hasOwnProperty("CP")){
        return criticalpairToPath(cr["CP"], handle_points, bp_scale, start_pos, unit_length);
    }else if (cr.hasOwnProperty("CR")){
        return confluentregionToPath(cr["CR"], handle_points, bp_scale, start_pos, unit_length);
    }else{
        console.log(`Error ${cr}`);
        return 1;
    }
};

const htgap = (read) => {
    let sum = 0;
    if (read.units[read.units.length-1].hasOwnProperty("G")){
        sum += read.units[read.units.length-1].G;
    }
    if (read.units[0].hasOwnProperty("G")){
        sum += read.units[0].G;
    }
    return sum;
};

const calcGap = (reads)=>{
    const len = reads.length;
    const sum = reads.map(read => htgap(read))
          .reduce((acc,x)=>acc+x, 0);
    return sum / len * 2;
};

const kFormatter = (num)=> {
    return Math.abs(num) > 999 ? Math.sign(num)*((Math.abs(num)/1000).toFixed(1)) + 'k' : Math.sign(num)*Math.abs(num);
};

const contigToHTML = (contig) =>{
    const start = kFormatter(contig["start_unit"]*150);
    const end = kFormatter(contig["end_unit"]*150);
    const direction = contig["direction"];
    return `<ul>
<li>Start:${start} bp</li>
<li>End:${end} bp</li>
<li>Direction:${direction} </li>
</ul>`;
};

const criticalpairToHTML = (cp,idx, reads) => {
    const count = reads.length;
    const meangap = calcGap(reads);
    const header = `<div>CriticalPair:${idx}</div>`;
    const contig1 = contigToHTML(cp["contig1"]);
    const contig2 = contigToHTML(cp["contig2"]);
    const support = `Supporing Reads:${count}<br>`;
    const gap = `Mean gap length:${meangap.toFixed(1)}`;
    return header + contig1 + contig2 + support + gap;
};

const confluentregionToHTML = (cr,idx, reads) => {
    const count = reads.length;
    const meangap = calcGap(reads);
    const header = `<div>ConfluentRegion:${idx}</div>`;
    const contig = contigToHTML(cr["pos"]);
    const support = `Supporing Reads:${count}<br>`;
    const gap = `Mean gap length:${meangap.toFixed(1)}`;
    return header + contig + support + gap;
};

const crToHTML = (cr,idx, reads) => {
    // Input: JSON object, Array
    // Output: String
    // Requirements: Critical region object
    // Return the HTML contents corresponds to the given cr.
    if (cr.hasOwnProperty("CP")){
        return criticalpairToHTML(cr["CP"],idx, reads);
    }else if (cr.hasOwnProperty("CR")){
        return confluentregionToHTML(cr["CR"],idx, reads);
    }else{
        console.log(`Error ${cr}`);
        return "Error";
    }
};

const calcCoverageOf = (reads, contigs)=>{
    // Input: List of JSON object, List of JSON object.
    // Output: List of JSON object
    // Requirements: An element of the first argument should be a JSON object having following
    // members: name => String, cluster => List of Integer, units => List of JSON Object.
    // Each unit is either {'G':Integer} or {'E':[Integer, Integer]}
    // An element of the second argument should be a JSON object having
    // name => String, id => Integer, length => integer, coverages => List of Integer,
    // start_stop => List of Integer
    // Specification: Each object in the output list should have the following elements:
    // id => integer
    // cov => list of integer
    let results = contigs.map(covs => {
        const len = covs.length;
        let coverage = Array.from({length:len}, (_) => 0);
        return {id: covs.id,
                length: covs.length,
                cov: coverage
               };});
    console.log("Start");
    for (const read of reads){
        for (const unit of read.units){
            if (unit.hasOwnProperty('E')){
                const c = unit.E[0];
                const s = unit.E[1];
                const e = unit.E[3];
                const start = Math.min(s,e);
                const end = Math.max(s,e);
                Array(end-start)
                    .fill(0)
                    .map((_,idx)=>idx+start)
                    .forEach((i)=> results[c].cov[i] +=1);
            }
        }
    }
    console.log("Done");
    return results;
};

const plotData = (dataset, repeats, unit_length) =>
      Promise.all([dataset, repeats]
                  .map(file => d3.json(file)))
      .then(([values, repeats]) => {
          // Unpack
          // This is array.
          const contigs = values.contigs;
          // This is also an array.
          // const reads = values.reads;
          // Or select reads as you like.
          const reads = values.reads;
          const critical_regions = values.critical_regions;
          const clusters = values.clusters;
          // Calculate coordinate.
          const bp_scale = calcScale(contigs);
          const coverage_scale = calcCovScale(contigs);
          const start_pos = calcStartPosition(contigs);
          const readnum_scale = calcReadNumScale(contigs);
          const handle_points = calcHandlePoints(start_pos);
          const contig_num = start_pos.length;
          const scales = {"bp_scale":bp_scale,
                          "coverage_scale":coverage_scale,
                          "start_pos": start_pos,
                          "readnum_scale":readnum_scale,
                          "handle_points":handle_points,
                          "start_pos": start_pos};
          // Draw contigs.
          contigs_layer
              .selectAll(".contig")
              .data(contigs)
              .enter()
              .append("path")
              .attr("class","contig")
              .attr("d", contig =>  {
                  const end = (contig.id == contig_num-1 ? Math.PI*2 : start_pos[contig.id+1]) - theta_margin;
                  const arc = d3.arc()
                        .innerRadius(contig_radius)
                        .outerRadius(contig_radius +  contig_thick)
                        .startAngle(start_pos[contig.id])
                        .endAngle(end);
                  return arc();
              })
              .attr("fill",c => d3.schemeCategory10[c.id% 10]);
          // Draw repeat.
          contigs_layer
              .selectAll(".repeats")
              .data(repeats.flatMap(rp => rp.reps))
              .enter()
              .append("path")
              .attr("class","repeats")
              .attr("d", repeat => {
                  const arc = d3.arc()
                        .innerRadius(contig_radius - 3)
                        .outerRadius(contig_radius + contig_thick + 3)
                        .startAngle(start_pos[repeat.id] + bp_scale(repeat.start))
                        .endAngle(start_pos[repeat.id] + bp_scale(repeat.end));
                  return arc();
              })
              .attr("fill", "gray");
          // Draw coverage
          coverage_layer
              .selectAll(".coverage")
              .data(contigs)
              .enter()
              .append("path")
              .attr("class","coverage")
              .attr("d", contig => {
                  const start = start_pos[contig.id];
                  const arc = d3.lineRadial()
                        .angle((_,i) => start + bp_scale(i * unit_length))
                        .radius(d => coverage_scale(d));
                  return arc(contig.coverages);
              })
              .attr("fill","none")
              .attr("stroke",c => d3.schemeCategory10[c.id% 10]);
          // Draw start/stop reads.
          start_stop_layer
              .selectAll(".start-stop-count")
              .data(contigs.flatMap(c => {
                  const start = start_pos[c.id] - Math.PI/2;
                  return c.start_stop.map((num,i) => {
                      const radian = start + bp_scale(i * unit_length);
                      const r = contig_radius + contig_thick/2;
                      const x = r * Math.cos(radian);
                      const y = r * Math.sin(radian);
                      return {"r":readnum_scale(num), "x": x, "y":y, "id":c.id};
                  });
              }))
              .enter()
              .append("circle")
              .attr("class",".start-stop-count")
              .attr("r", stst => stst.r)
              .attr("cx",stst => stst.x)
              .attr("cy",stst => stst.y)
              .attr("fill",stst => d3.schemeCategory10[stst.id % 10]);
          // Draw reads
          read_layer
              .selectAll(".read")
              .data(reads)
              .enter()
              .append("path")
              .attr("class","read")
              .attr("d",read => readToPath(read,handle_points,bp_scale,start_pos,unit_length))
              .attr("fill","none")
              .attr("opacity",0.3)
              .attr("stroke",read => "black");
          // Draw critical regions.
          cr_layer
              .selectAll(".cr")
              .data(critical_regions)
              .enter()
              .append("path")
              .attr("class", "cr")
              .attr("d", cr => crToPath(cr, handle_points, bp_scale, start_pos, unit_length))
              .attr("stroke", (cr,idx) => {
                  const searched = clusters.filter(cl => {
                      return cl.critical_regions.includes(idx);
                  })
                        .map((d,idx) => {
                            return {cluster:d,
                                    cluster_idx:idx};
                        });
                  if (searched.length == 0){
                      return "black";
                  }
                  const {cluster, cluster_idx} = searched[0];
                  return d3.schemeCategory10[(cluster_idx+1)%10];
              })
              .attr("stroke-width", (cr,idx) => (cr.hasOwnProperty("CP")) ? 5 : 100)
              .attr("stroke-linecap", (cr,idx) => (cr.hasOwnProperty("CP")) ? "round" : "none")
              .attr("opacity",cr => (cr.hasOwnProperty("CP")) ? 0.4 : 0.5)
              .attr("fill",  (cr,idx) => d3.schemeCategory10[(idx+1)%10])
              .on("mouseover", function(d,idx) {
                  // First, determine the cluster.
                  const searched = clusters.filter(cl => {
                      return cl.critical_regions.includes(idx);
                  })
                        .map((d,idx) => {
                            return {cluster:d,
                                    cluster_idx:idx};
                        });
                  if (searched.length == 0){
                      return;
                  }
                  const {cluster, cluster_idx} = searched[0];
                  const supporting_reads = reads.filter(r => r['cluster'] == cluster_idx);
                  tooltip.style("opacity", 0.9);
                  const contents = crToHTML(d,idx,supporting_reads);
                  tooltip.html(contents)
                      .style("left", (d3.event.pageX + 50) + "px")	
                      .style("top", (d3.event.pageY - 50) + "px");
                  temp_coverage_layer
                      .selectAll(".tempcoverage")
                      .data(cluster.contigs)
                      .enter()
                      .append("path")
                      .attr("class", "tempcoverage")
                      .attr("d", contig => {
                          console.log("Writing");
                          const start = start_pos[contig.id];
                          const arc = d3.lineRadial()
                                .angle((_,i) => start + bp_scale(i * unit_length))
                                .radius(d => coverage_scale(d));
                          return arc(contig.coverages);
                      })
                      .attr("fill","none")
                      .attr("opacity", 0.9)
                      .attr("stroke-width", 2)
                      .attr("stroke", d3.schemeCategory10[(idx +1)% 10]);
              })
              .on("mouseout", d => {
                  console.log("Detauch");
                  temp_coverage_layer
                      .selectAll(".tempcoverage")
                      .remove();
                  tooltip.style("opacity",0);
                  console.log("Detauched");
              });
          info.append("div")
              .attr("class","numofgapread")
              .append("p")
              .text(`Gap Read:${getNumOfGapRead(reads)} out of ${reads.length}`);
          // Draw ticks.
          const b_tick = svg.append("g")
                .attr("class","scale")
                .attr("transform",`translate(0,40)`);
          b_tick.append("text")
              .text("Base Pair Scale");
          {
              const bscale = d3.scaleLinear()
                    .domain([0,100000])
                    .range([contig_radius*bp_scale(0),contig_radius*bp_scale(100000)]);
              b_tick.append("g")
                  .attr("transform","translate(50,5)")
                  .call(d3.axisBottom(bscale)
                        .tickFormat(d3.format(".2s"))
                        .ticks(4));
          }
          const c_tick = svg.append("g")
                .attr("class","scale")
                .attr("transform",`translate(0,100)`);
          c_tick.append("text")
              .text("Coverage Scale");
          {
              const cscale = d3.scaleLinear()
                    .domain([0,1500])
                    .range([0,coverage_scale(1500)-coverage_scale(0)]);
              c_tick.append("g")
                  .attr("transform",`translate(50,5)`)
                  .call(d3.axisBottom(cscale)
                        .tickFormat(d3.format(".2s"))
                        .ticks(4));
          }
          const g_tick = svg.append("g")
                .attr("class","scale")
                .attr("transform",`translate(0,160)`);
          g_tick.append("text")
              .text("Gap Scale");
          {
              const gscale = d3.scaleLog()
                    .domain([gap_min,5*gap_max])
                    .range([0, 5*(gap_max_radius-gap_min_radius)]);
              g_tick.append("g")
                  .attr("transform",`translate(50,5)`)
                  .call(d3.axisBottom(gscale)
                        .tickFormat(d3.format(".2s"))
                        .ticks(1)
                       );
          }
          const n_tick = svg.append("g")
                .attr("class","scale")
                .attr("transform", `translate(0,220)`);
          n_tick.append("text")
              .text("Number of Reads");
          {
              const sizes = [3,9,20];
              n_tick.append("g")
                  .attr("transform",`translate(60,15)`)
                  .selectAll("specimen")
                  .data(sizes)
                  .enter()
                  .append("circle")
                  .attr("class","specimen")
                  .attr("cx", (_,i) => 20 *i)
                  .attr("cy", 0)
                  .attr("r" , r => readnum_scale(r))
                  .attr("fill","black");
              n_tick.append("g")
                  .attr("transform",`translate(60,35)`)
                  .selectAll("ticks")
                  .data(sizes)
                  .enter()
                  .append("text")
                  .attr("x", (_,i) => 20 *i)
                  .attr("y", 0)
                  .text(r => r);
          }
          return scales;
      })
      .then(ok => ok,
            why => console.log(why));

