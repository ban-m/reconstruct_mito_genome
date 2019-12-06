const width = 1000;
const height = 1000;
// Margin in radian
const theta_margin = 0.15;
const gap_position = 0.05;
// the height of a contigs.
const contig_thick = 10; 
const coverage_thick = 5;
const gap_jitters = d3.randomNormal(0,0.01);
const read_thick = 4;
const eplison = 5.001;
const jitters = d3.randomNormal(0,eplison);

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
const start_stop_layer = svg.append("g")
      .attr("transform", `translate(${width/2},${height/2+offset})`)
      .attr("class","start-stop-read");
const read_layer = svg.append("g")
      .attr("transform", `translate(${width/2},${height/2+offset})`)
      .attr("class","read");
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
    // should have either "Gap" or "Encode"
    let path = d3.path();
    let units = Array.from(read.units).reverse();
    const r = read_radius; // + jitters();
    let gap = 0;
    let unit = {};
    while(!unit.hasOwnProperty("Encode")){
        unit = units.pop();
        if (unit == undefined){
            return "";
        }else if (unit.hasOwnProperty("Gap")){
            gap = unit.Gap;
        }
    };
    // Current ID of the contig 
    let contig = unit.Encode[0];
    let start = start_pos[contig] - Math.PI/2;
    let radian = start + bp_scale(unit_length*unit.Encode[1]);
    if (gap != 0){
        path.moveTo(gap_scale(gap) * Math.cos(radian), gap_scale(gap)*Math.sin(radian));
        path.lineTo(read_radius * Math.cos(radian), read_radius * Math.sin(radian));
    }else{
        path.moveTo(read_radius * Math.cos(radian), read_radius * Math.sin(radian));
    }
    // gap = 0;
    for (unit of units.reverse()){
        if (unit.hasOwnProperty("Gap")){
            // if (unit.Gap > unit_length * 2){
            //     gap = unit.Gap;
            // }
            continue;
        }
        if (unit.Encode[0] == contig){
            radian = start + bp_scale(unit_length*unit.Encode[1]);
            path.lineTo(r * Math.cos(radian), r*Math.sin(radian));
        }else{
            // Change contig. Connect them.
            // If there are remaining gap, clean them.
            // if (gap != 0){
            //     const control_radian = start_pos[contig] - Math.PI/2;
            //     const new_radian = control_radian - gap_position;
            //     const control_x = handle_points_radius * Math.cos(control_radian);
            //     const control_y = handle_points_radius * Math.sin(control_radian);
            //     const jt = gap_jitters();
            //     path.quadraticCurveTo(control_x, control_y, r * Math.cos(new_radian), r * Math.sin(new_radian));
            //     path.moveTo(gap_scale(gap) * Math.cos(new_radian + jt), gap_scale(gap)*Math.sin(new_radian + jt));
            //     path.lineTo(r * Math.cos(new_radian), r * Math.sin(new_radian));
            // }
            // gap = 0;
            const new_radian = start_pos[unit.Encode[0]];
            radian = new_radian + bp_scale(unit_length*unit.Encode[1]) - Math.PI/2;
            // Bezier Curve to new point from here.
            const control_radius = handle_points[contig][unit.Encode[0]] - Math.PI/2;
            const control_x = handle_points_radius*Math.cos(control_radius);
            const control_y = handle_points_radius*Math.sin(control_radius);
            contig = unit.Encode[0];
            start = start_pos[contig] - Math.PI/2;
            path.quadraticCurveTo(control_x,control_y,r*Math.cos(radian),r*Math.sin(radian));
        }
    }
    return path.toString();
};

const calcID = (read,unit_length)=>{
    // Input: Json object
    // Output: JSON object having "type" property and "id" property(maybe).
    // Requirements: read should have "units" property, which is a vector
    // and each of element should have eigther "Gap" or "Encode" type.
    // Returns the most assigned type of given read.
    const gap = read
          .units
          .filter(unit => unit.hasOwnProperty("Gap"))
          .reduce((g, unit) => g + unit.Gap,0);
    const summary = read
          .units
          .filter(unit => unit.hasOwnProperty("Encode"))
          .map(unit => unit.Encode[0])
          .reduce((map,ctg)=>{
              if (map.has(ctg)){
                  map.set(ctg,map.get(ctg)+unit_length);
              }else{
                  map.set(ctg,unit_length);
              }
              return map;
          }, new Map());
    let max = undefined;
    summary
        .forEach((len,ctg)=>{
            if (max == undefined || max.len < len){
                max = {"ctg":ctg, "len":len};
            }else {
                max = max;
            }});
    if (max == undefined){
        return {"type":"Gap"};
    }else{
        return (max.len < gap ? {"type":"Gap"} : {"type":"Contig", "id":max.ctg});
    }
};


const selectRead = (read,unitlen) => {
    // Input: JSON object, Num
    // Output: boolean
    // Requirements: input object should have "units" property,
    // which is actually vector of object with "Gap" or "Encode" property.
    // Filter read as you like.
    const from = 0;
    const to = 1;
    const set = new Set(read.units.filter(u => u.hasOwnProperty("Encode")).map(u => u.Encode[0]));
    const max_gap = Math.max(...read.units.filter(u => u.hasOwnProperty("Gap")).map(u => u.Gap));
    return true;
    // return set.has(4) && set.size == 1;
    // return set.has(from) && set.has(to) && read.units.length > 15 ;
    // return read.units.length < 140 && read.units.length > 75 && set.size > 1 && set.has(0) && set.has(1) && max_gap < 4000;
    // return set.size == 1 && (set.has(0) || set.has(1)) && calcID(read,unitlen).type == "Contig" ; // && max_gap < 4000;
    // return set.has(0) && set.has(4) ;
    // return !set.has(6);
};

const getNumOfGapRead = reads => {
    // Input: [JSON object]
    // Output: Num
    // Requirements: each element should be 'read' object.
    // Return numbers of reads which is just Gap.
    return reads.filter(read => {
        let units = Array.from(read.units);
        let unit = {};
        while(!unit.hasOwnProperty("Encode")){
            unit = units.pop();
            if (unit == undefined){
                return true;
            }
        };
        return false;
    }).length;
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
          const reads = values.reads.filter(r => selectRead(r,unit_length));
          // let reads = values.reads.slice(0,10);
          // reads.push({"name":"test",
          //             "units":[{"Gap":1000},
          //                      {"Encode":[0,0]},{"Encode":[0,1]},{"Encode":[0,2]},{"Encode":[2,100]},
          //                      {"Gap":2000}]});
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
              .attr("stroke",read => {
                  const identity = calcID(read,unit_length);
                  if (identity.type == "Gap"){
                      return "black";
                  }else{
                      return d3.schemeCategory10[identity.id % 10];
                  }
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
          // const pic = document.getElementById("plot");
          //get svg source.
          // var serializer = new XMLSerializer();
          // var source = serializer.serializeToString(pic);
          // var svgBlob = new Blob([source], {type:"image/svg+xml;charset=utf-8"});
          // var svgUrl = URL.createObjectURL(svgBlob);
          // var downloadLink = document.createElement("a");
          // downloadLink.href = svgUrl;
          // downloadLink.download = "newesttree.svg";
          // document.body.appendChild(downloadLink);
          // downloadLink.click();
          // document.body.removeChild(downloadLink);
          return scales;
      })
      .then(ok => ok,
            why => console.log(why));

const make_path_between = (cr, scales,unit_length)=>{
    // Input: JSON object, JSON object
    // Output: String
    // Requirements: Critical region object, scales
    // Return the path btw critical region, or just the point.
    let p = d3.path();
    const inner = cr.CP;
    const start1 = scales.start_pos[inner.contig1.contig] +
          scales.bp_scale(inner.contig1.start_unit * unit_length) - Math.PI/2;
    const end1 = scales.start_pos[inner.contig1.contig] +
          scales.bp_scale(inner.contig1.start_unit * unit_length) - Math.PI/2;
    const start2 = scales.start_pos[inner.contig2.contig] +
          scales.bp_scale(inner.contig2.start_unit * unit_length) - Math.PI/2;
    const end2 = scales.start_pos[inner.contig2.contig] +
          scales.bp_scale(inner.contig2.start_unit * unit_length) - Math.PI/2;
    const hp = scales.handle_points[inner.contig1.contig][inner.contig2.contig] - Math.PI/2;
    const hpx = handle_points_radius * Math.cos(hp);
    const hpy = handle_points_radius * Math.sin(hp);
    p.moveTo(read_radius * Math.cos(start1), read_radius * Math.sin(start1));
    p.lineTo(read_radius * Math.cos(start2), read_radius * Math.sin(start2));
    // p.arc(0,0,read_radius, start1, end1);
    // p.quadraticCurveTo(hpx,hpy,read_radius*Math.cos(start2), read_radius*Math.sin(start2));
    // p.arc(0,0,read_radius, start2, end2);
    // p.quadraticCurveTo(hpx,hpy, read_radius*Math.cos(start1), read_radius*Math.sin(start1));
    // p.closePath();
    return p.toString();
};

const overlay_cr = (scales,critical_regions, unit_length) =>
      d3.json(critical_regions)
      .then(critical_regions => {
          console.log(scales);
          console.log(critical_regions);
          contigs_layer
              .selectAll("critical_region")
              .data(critical_regions.filter(d => d.hasOwnProperty("CP")))
              .enter()
              .append("path")
              .attr("class","critical_region")
              .attr("d", cr => make_path_between(cr, scales,unit_length))
              .attr("stroke-width",4)
              .attr("stroke", "yellow");
          contigs_layer
              .selectAll("critical_region")
              .data(critical_regions.filter(d => d.hasOwnProperty("RJ")))
              .enter()
              .append("path")
              .attr("class","critical_region")
              .attr("d", d => {
                  const inner = d.RJ.pos;
                  const r = scales.start_pos[inner.contig];
                  const start = r+ scales.bp_scale(unit_length*inner.start_unit);
                  const end = r + scales.bp_scale(unit_length*inner.end_unit);
                  const arc = d3.arc()
                        .innerRadius(read_radius - 10)
                        .outerRadius(read_radius)
                        .startAngle(start)
                        .endAngle(end);
                  return arc();
              })
              .attr("fill", "yellow");
      })
      .then(ok => console.log("OK"),
            why => console.log(`Error:${why}`));
