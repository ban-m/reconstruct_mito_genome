const width = 1000;
const height = 1000;
// Margin in radian
const theta_margin = 0.10;
// the height of a contigs.
const contig_thick = 10; 
const coverage_thick = 5;
const read_thick = 4;
const eplison = 0.7;

// Radius
const contig_radius = 350;
const coverage_min = contig_radius+contig_thick;
const coverage_max = 450;
const handle_points_radius = 100;
const read_radius = contig_radius-30;
const gap_min_radius = read_radius;
const gap_max_radius = contig_radius-3;
const gap_min = 100;
const gap_max = 5000;
const gap_scale = d3.scaleLog()
      .domain([gap_min,gap_max])
      .range([gap_min_radius, gap_max_radius])
      .clamp(true);

const svg = d3.select("#plot")
      .append("svg")
      .attr("width",width)
      .attr("height",height);
const contigs_layer = svg.append("g")
      .attr("transform", `translate(${width/2},${height/2})`)
      .attr("class","contigs");
const coverage_layer = svg.append("g")
      .attr("transform", `translate(${width/2},${height/2})`)
      .attr("class","coverages");
const read_layer = svg.append("g")
      .attr("transform", `translate(${width/2},${height/2})`)
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
    const min = Math.min(...contigs.map(contig => Math.min(...contig.coverages)));
    return d3.scaleLinear()
        .domain([min,max])
        .range([coverage_min,coverage_max]);
};

const readToPath = (read,handle_points,bp_scale,start_pos,unit_length)=>{
    // Input: JSON object, Array[Array[Num]], d3.scale, Array[Num], Num
    // Output: String
    // Requirements: read should have units attribute, each of which elements
    // should have either "Gap" or "Encode"
    let path = d3.path();
    let units = Array.from(read.units).reverse();
    let jitters = d3.randomNormal(0,eplison);
    let gap = 0;
    let unit = {};
    while(!unit.hasOwnProperty("Encode")){
        unit = units.pop();
        if (unit == undefined){
            console.log(`${read} is full of gap!`);
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
    while (1){
        let next = units.pop();
        if (next == undefined){
            break;
        }else if (next.hasOwnProperty("Gap")){
            continue;
        };
        if (next.Encode[0] == contig){
            radian = start + bp_scale(unit_length*next.Encode[1]);
            const r = read_radius + jitters();
            path.lineTo(r * Math.cos(radian), r*Math.sin(radian));
        }else{
            // Change contig. Connect them.
            const new_radian = start_pos[next.Encode[0]];
            radian = new_radian + bp_scale(unit_length*next.Encode[1]) - Math.PI/2;
            // Bezier Curve to new point from here.
            const control_radius = handle_points[contig][next.Encode[0]] - Math.PI/2;
            const control_x = handle_points_radius*Math.cos(control_radius);
            const control_y = handle_points_radius*Math.sin(control_radius);
            contig = next.Encode[0];
            start = start_pos[contig] - Math.PI/2;
            const r = read_radius + jitters();
            path.quadraticCurveTo(control_x,control_y,r*Math.cos(radian),r*Math.sin(radian));
        }
        unit = next;
    }
    return path.toString();
};

const calcID = (read,unit_length)=>{
    // Input: Json object
    // Output: Num
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


const selectRead = read => {
    // Input: JSON object
    // Output: boolean
    // Requirements: input object should have "units" property,
    // which is actually vector of object with "Gap" or "Encode" property.
    // Filter read as you like.
    const from = 0;
    const to = 1;
    const set = new Set(read.units.filter(u => u.hasOwnProperty("Encode")).map(u => u.Encode[0]));
    return true;
    // return set.has(from) && set.has(to) && read.units.length > 15 ;
    // return set.has(2) &&  set.has(1) && set.has(0) && read.units.length > 140;
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
          const reads = values.reads.filter(selectRead);
          // let reads = values.reads.slice(0,10);
          // reads.push({"name":"test",
          //             "units":[{"Gap":1000},
          //                      {"Encode":[0,0]},{"Encode":[0,1]},{"Encode":[0,2]},{"Encode":[2,100]},
          //                      {"Gap":2000}]});
          // Calculate coordinate.
          const bp_scale = calcScale(contigs);
          const coverage_scale = calcCovScale(contigs);
          const start_pos = calcStartPosition(contigs);
          const handle_points = calcHandlePoints(start_pos);
          const contig_num = start_pos.length;
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
          // console.log(reads);
          // Draw reads
          read_layer
              .selectAll(".read")
              .data(reads)
              .enter()
              .append("path")
              .attr("class","read")
              .attr("d",read => readToPath(read,handle_points,bp_scale,start_pos,unit_length))
              .attr("fill","none")
              .attr("opacity",0.2)
              .attr("stroke",read => {
                  const identity = calcID(read,unit_length);
                  if (identity.type == "Gap"){
                      return "black";
                  }else{
                      return d3.schemeCategory10[identity.id % 10];
                  }
              });
          // Draw ticks.
          const c_tick = svg.append("g")
                .attr("class","scale")
                .attr("transform",`translate(0,100)`);
          c_tick.append("text")
              .text("Coverage Scale");
          c_tick.append("g")
              .attr("transform",`translate(-200,0)`)
              .call(d3.axisBottom(coverage_scale)
                    .ticks(2));
          const g_tick = svg.append("g")
                .attr("class","scale")
                .attr("transform",`translate(0,150)`);
          g_tick.append("text")
              .text("Gap Scale");
          g_tick.append("g")
              .attr("transform",`translate(-200,0)`)
              .call(d3.axisBottom(gap_scale)
                    .ticks(1));
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
      })
      .then(ok => console.log("OK"),
            why => console.log(why));
