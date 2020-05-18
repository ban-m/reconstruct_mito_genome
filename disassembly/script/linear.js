const left_contig_margin = 200;
const total_reference_length = 600;
const reference_margin = 10;
const reference_label_margin = 30;
const contig_margin = 70;
const contig_thick = 10;
const contig_label_margin = 20;
const coverage_thick = 5;
const coverage_interval = 200;
const read_thick = 4;
const eplison = 5.001;
const confluent_margin = 0.01;
const coverage_min = contig_thick;
const coverage_max = 100;
const handle_points_radius = 100;
// Circle radius. Used in start/stop circle.
const min_radius = 1;
const max_radius = 8;
const min_read_num = 2;
// Gaps at the head/tail of the read.
const gap_min = 2;
const gap_max = 20;
const gap_min_size = 500;
const gap_max_size = 2000;
const width = 1000;
const tick_top_margin = 20;
const tick_margin = 50;
const tick_interval = 100000;
const top_margin = 50;
const master_circle = 20;

const gap_scale = d3
  .scaleLog()
  .domain([gap_min_size, gap_max_size])
  .range([gap_min, gap_max])
  .clamp(true);

const calcScale = (contigs) => {
  // Input: Array of JSON object
  // Output: d3 Scale object
  // Requirements: each input object should have "length" attribute
  const sum = contigs.map((c) => c.length).reduce((x, y) => x + y);
  return d3.scaleLinear().domain([0, sum]).range([0, total_reference_length]);
};

const calcStartPosition = (contigs, bp_scale) => {
  // Input: Array of JSON object, d3 scale object
  // Output: Array[Num]
  // Requirements: each input object should have "length" attribute
  // Map from the index(id) of contig into the start position of the contig(the y-axis).
  // Note that the start position may be negative, as the center is (left_margin, height/2)
  let start_pos = [];
  let cum = 0;
  contigs.forEach((contig, idx) => {
    start_pos.push(cum);
    cum += bp_scale(contig.length) + contig_margin;
  });
  return start_pos;
};

const calcStartMC = (references, bp_scale) => {
  let pos = new Map();
  let cum = 0;
  for (const reference of references) {
    const prev = cum;
    cum += bp_scale(reference.length) + reference_margin;
    pos.set(reference.id, [prev, cum]);
  }
  return pos;
};

const calcHandlePoints = (start_pos) => {
  // Input: Array[Num]
  // Output: Array[[Num => Num]]
  // Requirement: None
  // Map from combinations of ID to the handle points of them.
  let handle_points = new Array();
  for (let i = 0; i < start_pos.length; i++) {
    handle_points.push(new Array(start_pos.length));
  }
  start_pos.forEach((v1, k1) => {
    start_pos.forEach((v2, k2) => {
      if (k1 == k2) {
        handle_points[k1][k2] = v1 + contig_margin / 2;
      } else {
        handle_points[k1][k2] = (v1 + v2) / 2;
      }
    });
  });
  return handle_points;
};

const calcCovScale = (contigs) => {
  // Input: Array on JSON object
  // Output: d3.scale object
  // Requirements: each input object should have "length" attribute
  // Scale for convert coverage into radius.
  const max = Math.max(
    ...contigs.map((contig) => Math.max(...contig.coverages))
  );
  // const min = Math.min(...contigs.map(contig => Math.min(...contig.coverages)));
  return d3.scaleLinear().domain([0, max]).range([coverage_min, coverage_max]);
};

const calcReadNumScale = (contigs) => {
  // Input: Array on JSON object
  // Output: d3.scale object
  // Requirements: Each object in the argument should have an array of integer, which is
  // named "start_stop."
  // Calculate the scale for start/stop vizualization.
  const total = contigs.flatMap((c) => c.start_stop).reduce((x, y) => x + y);
  const num = contigs.map((c) => c.start_stop.length).reduce((x, y) => x + y);
  const max = Math.max(...contigs.flatMap((c) => c.start_stop));
  return d3
    .scaleLog()
    .domain([min_read_num, max])
    .range([min_radius, max_radius])
    .clamp(true);
};

const readToPath = (read, handle_points, bp_scale, start_pos, unit_length) => {
  // Input: JSON object, Array[Array[Num]], d3.scale, Array[Num], Num
  // Output: String
  // Requirements: read should have units attribute, each of which elements
  // should have either "G"(for Gap) or "E"(for Encode)
  let path = d3.path();
  let gap = 0;
  let index = 0;
  while (read.is_unit[index] == 0 && index < read.is_unit.length) {
    gap = read.contig[index];
    index += 1;
  }
  // Current ID of the contig
  let contig = read.contig[index];
  let current_unit = read.unit[index];
  let y = start_pos[contig] + bp_scale(unit_length * current_unit);
  let x = 0;
  if (gap != 0) {
    path.moveTo(x - gap_scale(gap), y);
    path.lineTo(x, y);
  } else {
    path.moveTo(x, y);
  }
  for (; index < read.is_unit.length; index += 1) {
    if (read.is_unit[index] == 0) {
      continue;
    }
    const diff = Math.abs(read.unit[index] - current_unit);
    if (read.unit[index] == contig && diff < 50) {
    } else {
      y = start_pos[contig] + bp_scale(unit_length * current_unit);
      path.lineTo(x, y);
      // Change contig. Connect them.
      const new_y =
        start_pos[read.unit[index]] + bp_scale(unit_length * read.unit[index]);
      const new_x = 0;
      // Bezier Curve to new point from here.
      const control_y = (y + new_y) / 2;
      const control_x = 10; //handle_points[contig][unit.E[0]];
      contig = read.contig[index];
      path.quadraticCurveTo(control_x, control_y, new_x, new_y);
      x = new_x;
      y = new_y;
    }
    current_unit = read.unit[index];
  }
  return path.toString();
};

// Input: [JSON object]
// Output: Num
// Requirements: each element should be 'read' object.
// Return numbers of reads which is just Gap.
const getNumOfGapRead = (reads) =>
  reads.filter((read) => read.is_unit.every((x) => x == 0)).length;

const htgap = (read) => {
  let sum = 0;
  if (read.is_unit[read.is_unit.length - 1] == 0) {
    sum += read.contig[read.contig.length - 1];
  }
  if (read.is_unit[0] == 0) {
    sum += read.contig[0];
  }
  return sum;
};

const calcGap = (reads) =>
  (reads.map(htgap).reduce((acc, x) => acc + x) / reads.length) * 2;

const kFormatter = (num) => {
  return Math.abs(num) > 999
    ? Math.sign(num) * (Math.abs(num) / 1000).toFixed(1) + "k"
    : Math.sign(num) * Math.abs(num);
};

const contigToHTML = (contig) => {
  const start = kFormatter(contig["start_unit"] * 150);
  const end = kFormatter(contig["end_unit"] * 150);
  const direction = contig["direction"];
  return `<ul>
<li>Start:${start} bp</li>
<li>End:${end} bp</li>
<li>Direction:${direction} </li>
</ul>`;
};

const compare_name = (contig1, contig2) => {
  const name1_vs = contig1.name.split("_").map((c) => parseInt(c));
  const name2_vs = contig2.name.split("_").map((c) => parseInt(c));
  if (name1_vs[0] < name2_vs[0]) {
    return -1;
  } else if (name1_vs[0] > name2_vs[0]) {
    return +1;
  } else {
    return name1_vs[1] - name2_vs[1];
  }
};

const convertRead = (read) => {
  // Input:Read
  // Output: An object consists of three typed array: one Uint8Array and two Uint16Array.
  // They are named "is_unit", "contig", and "unit" respectively.
  // 1 is true, 0 is false for the is_unit array.
  const name = read.name.repeat(1);
  const cluster = read.cluster;
  const is_unit = Uint8Array.from(
    read.units.map((u) => (u.hasOwnProperty("E") ? 1 : 0))
  );
  const contig = Uint16Array.from(
    read.units.map((u) => (u.hasOwnProperty("E") ? u.E[0] : u.G))
  );
  const unit = Uint16Array.from(
    read.units.map((u) => (u.hasOwnProperty("E") ? u.E[1] : 0))
  );
  return {
    name: name,
    cluster: cluster,
    is_unit: is_unit,
    contig: contig,
    unit: unit,
  };
};

const plotData = (dataset, alignments, unit_length) =>
  Promise.all([dataset, alignments].map((file) => d3.json(file)))
    .then(([values, contig_alns]) => {
      // Margin between contigs
      const contigs = values.contigs; //.sort(compare_name);
      const references = contig_alns.references;
      const colorScheme = d3.schemePaired;
      const bp_scale = calcScale(references);
      const height =
        contigs
          .map((c) => bp_scale(c.length) + contig_margin)
          .reduce((x, y) => x + y) + top_margin;
      const svg = d3
        .select("#plot")
        .append("svg")
        .attr("width", width)
        .attr("height", height);
      const contigs_layer = svg
        .append("g")
        .attr("transform", `translate(${left_contig_margin},${top_margin})`)
        .attr("class", "contigs");
      const coverage_layer = svg
        .append("g")
        .attr("transform", `translate(${left_contig_margin},${top_margin})`)
        .attr("class", "coverages");
      const start_stop_layer = svg
        .append("g")
        .attr("transform", `translate(${left_contig_margin},${top_margin})`)
        .attr("class", "start-stop-read");
      const read_layer = svg
        .append("g")
        .attr("transform", `translate(${left_contig_margin},${top_margin})`)
        .attr("class", "read");
      const aln_layer = svg
        .append("g")
        .attr("transform", `translate(${left_contig_margin},${top_margin})`)
        .attr("class", "critical-region");
      const info = d3.select("#info");
      // This is also an array.
      // const reads = values.reads;
      // Or select reads as you like.
      const reads = values.reads.map(convertRead);
      // const critical_regions = [values.critical_regions[selected_region]];
      const clusters = values.clusters;
      // Calculate coordinate.
      const coverage_scale = calcCovScale(contigs);
      const start_pos = calcStartPosition(contigs, bp_scale);
      const readnum_scale = calcReadNumScale(contigs);
      const handle_points = calcHandlePoints(start_pos);
      const contig_num = start_pos.length;
      const scales = {
        bp_scale: bp_scale,
        coverage_scale: coverage_scale,
        start_pos: start_pos,
        readnum_scale: readnum_scale,
        handle_points: handle_points,
        start_pos: start_pos,
      };
      // Draw coverage
      const each_coverage_layer = coverage_layer
        .selectAll(".coverage")
        .data(contigs)
        .enter()
        .append("g")
        .attr("class", (c) => `coverage-${c.id}`)
        .attr("transform", (contig) => `translate(0,${start_pos[contig.id]})`);
      each_coverage_layer
        .append("path")
        .attr("class", "coverage")
        .attr("d", (contig) => {
          const path = d3
            .line()
            .x((cov) => -coverage_scale(cov))
            .y((_, i) => bp_scale(i * unit_length));
          return path(contig.coverages);
        })
        .attr("fill", "none")
        .attr("stroke", (c) => d3.schemeCategory10[c.id % 10]);
      each_coverage_layer
        .each((c) => {
          const max_cov = Math.max(coverage_interval, Math.max(...c.coverages));
          const domain = [0, max_cov];
          const range = [0, -coverage_scale(domain[1])];
          const local_bp_scale = d3.scaleLinear().domain(domain).range(range);
          const tV = Array.from({
            length: Math.floor(max_cov / coverage_interval),
          }).map((_, idx) => (idx + 1) * coverage_interval);
          d3.selectAll(`.coverage-${c.id}`).call(
            d3
              .axisTop(local_bp_scale)
              .tickFormat(d3.format(".2s"))
              .tickValues(tV)
          );
        })
        .selectAll(".tick text")
        .attr("transform", "rotate(-15)");
      // Draw reads
      // read_layer
      //   .selectAll(".read")
      //   .data(reads)
      //   .enter()
      //   .append("path")
      //   .attr("class", "read")
      //   .attr("d", (read) =>
      //         readToPath(read, handle_points, bp_scale, start_pos, unit_length)
      //   )
      //   .attr("fill", "none")
      //   .attr("opacity", 0.3)
      //   .attr("stroke", (read) => "black");
      info
        .append("div")
        .attr("class", "numofgapread")
        .append("p")
        .text(`Gap Read:${getNumOfGapRead(reads)} out of ${reads.length}`);
      // Draw contigs
      const selection_of_each_contig = contigs_layer
        .selectAll(".contig")
        .data(contigs)
        .enter()
        .append("g")
        .attr("class", "contig")
        .attr("transform", (c) => `translate(0,${start_pos[c.id]})`);
      selection_of_each_contig
        .append("g")
        .attr("class", (c) => `contig-${c.id}`)
        .each((c) => {
          const domain = [0, c.length];
          const range = [0, bp_scale(c.length)];
          const local_bp_scale = d3.scaleLinear().domain(domain).range(range);
          const tV = Array.from({
            length: Math.floor(c.length / tick_interval),
          }).map((_, idx) => (idx + 1) * tick_interval);
          d3.selectAll(`.contig-${c.id}`).call(
            d3
              .axisLeft(local_bp_scale)
              .tickFormat(d3.format(".2s"))
              .tickValues(tV)
          );
        })
        .append("text")
        .attr(
          "transform",
          (c) =>
            `translate(-${coverage_max}, ${bp_scale(c.length / 2)}) rotate(90)`
        )
        .attr("color", "black")
        .attr("fill", "black")
        .text((c) => `${c.id}`);
      const reference_start = calcStartMC(references, bp_scale);
      selection_of_each_contig
        .selectAll(".mc")
        .data(references)
        .enter()
        .append("g")
        .attr(
          "transform",
          (m) => `translate(${reference_start.get(m.id)[0]},0)`
        )
        .attr("class", (_, idx) => `scale-${idx}`)
        .each((m, idx) => {
          const domain = [0, m.length];
          const range = [0, bp_scale(m.length)];
          const local_bp_scale = d3.scaleLinear().domain(domain).range(range);
          const tV = Array.from({
            length: Math.floor(m.length / tick_interval),
          }).map((_, idx) => (idx + 1) * tick_interval);
          d3.selectAll(`.scale-${idx}`).call(
            d3
              .axisTop(local_bp_scale)
              .tickFormat(d3.format(".2s"))
              .tickValues(tV)
          );
        })
        .append("text")
        .attr("x", (m) => bp_scale(m.length) / 2)
        .attr("y", -reference_label_margin)
        .attr("fill", "black")
        .text((m) => `${m.id}`);
      // Draw alignments
      aln_layer
        .selectAll(".alignments")
        .data(contig_alns.alignments)
        .enter()
        .append("line")
        .attr("x1", (aln) => {
          const ref_name = aln.reference_name;
          const start = bp_scale(Math.min(aln.ref_start, aln.ref_end));
          return reference_start.get(ref_name)[0] + start;
        })
        .attr("x2", (aln) => {
          const ref_name = aln.reference_name;
          const end = bp_scale(Math.max(aln.ref_start, aln.ref_end));
          return reference_start.get(ref_name)[0] + end;
        })
        .attr("y1", (aln) => {
          const elm = contigs.find((c) => c.name === aln.query_name);
          const start = bp_scale(Math.min(aln.query_start, aln.query_end));
          return start_pos[elm.id] + start;
        })
        .attr("y2", (aln) => {
          const elm = contigs.find((c) => c.name === aln.query_name);
          const end = bp_scale(Math.max(aln.query_start, aln.query_end));
          return start_pos[elm.id] + end;
        })
        .attr("stroke", "black")
        .attr("stroke_width", 2);
      return scales;
    })
    .then(
      (ok) => ok,
      (why) => console.log(why)
    );
