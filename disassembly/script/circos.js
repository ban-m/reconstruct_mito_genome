const width = 1000;
const height = 1000;
// Margin in radian
const theta_margin = 0.1;
// the height of a contigs.
const reference_tick = 10000;
const reference_break = 50000;
const gap_jitters = d3.randomNormal(0, 0.01);
const read_thick = 4;
const eplison = 5.001;
const jitters = d3.randomNormal(0, eplison);
const confluent_margin = 0.01;
// Radius
const contig_radius = 400;
const contig_thick = 10;
const coverage_min = contig_radius;
const coverage_max = coverage_min + 100;
const handle_points_radius = 100;
const read_radius = contig_radius - contig_thick;
const start_stop_radius = read_radius;
const start_stop_max = read_radius - 80;
const start_stop_break = 50;
const gap_thr = 2000;

// Circle radius
const offset = 0;
const svg = d3
  .select("#plot")
  .append("svg")
  .attr("width", width)
  .attr("height", height);
const contigs_layer = svg
  .append("g")
  .attr("transform", `translate(${width / 2},${height / 2 + offset})`)
  .attr("class", "contigs");
const coverage_layer = svg
  .append("g")
  .attr("transform", `translate(${width / 2},${height / 2 + offset})`)
  .attr("class", "coverages");
const temp_coverage_layer = svg
  .append("g")
  .attr("transform", `translate(${width / 2},${height / 2 + offset})`)
  .attr("class", "temp-coverages");
const start_stop_layer = svg
  .append("g")
  .attr("transform", `translate(${width / 2},${height / 2 + offset})`)
  .attr("class", "start-stop");
const read_layer = svg
  .append("g")
  .attr("transform", `translate(${width / 2},${height / 2 + offset})`)
  .attr("class", "read");
const cr_layer = svg
  .append("g")
  .attr("transform", `translate(${width / 2},${height / 2 + offset})`)
  .attr("class", "critical-region");
const tooltip = d3
  .select("body")
  .append("div")
  .attr("class", "tooltip")
  .style("opacity", 0);
const info = d3.select("#info");
const cr_info = d3.select("#cr-info");
const calcScale = (contigs) => {
  // Input: Array of JSON object
  // Output: d3 Scale object
  // Requirements: each input object should have "length" attribute
  // Convert base pair into radian
  const num = contigs.length;
  const total = contigs.map((c) => c.length).reduce((x, y) => x + y);
  return d3
    .scaleLinear()
    .domain([0, total])
    .range([0, 2 * Math.PI - num * theta_margin]);
};

const calcStartPosition = (contigs) => {
  // Input: Array of JSON object.
  // Output: Array[Num]
  // Requirements: each input object should have "length" attribute
  // Map from the index(id) of contig into the start position of the contig(in radian).
  const scale = calcScale(contigs);
  const max = Math.max(...contigs.map((c) => c.id));
  let cum_sum = 0;
  let start_pos = new Array(max);
  for (const contig of contigs) {
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
  for (let i = 0; i < start_pos.length; i++) {
    handle_points.push(new Array(start_pos.length));
  }
  const max = start_pos.length - 1;
  start_pos.forEach((v1, k1) => {
    start_pos.forEach((v2, k2) => {
      const next1 = k1 == max ? Math.PI * 2 - theta_margin : start_pos[k1];
      const next2 = k2 == max ? Math.PI * 2 - theta_margin : start_pos[k2];
      const val = (next1 + next2 + v1 + v2) / 4 - theta_margin / 2;
      handle_points[k1][k2] = val;
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
    .scaleLinear()
    .domain([0, max])
    .range([start_stop_radius, start_stop_max])
    .clamp(true);
};

const readToPath = (
  read,
  handle_points,
  bp_scale,
  start_pos,
  unit_length,
  ends
) => {
  // Input: JSON object, Array[Array[Num]], d3.scale, Array[Num], Num, Array[Num]
  // Output: String
  // Requirements: read should have units attribute, each of which elements
  // should have either "G"(for Gap) or "E"(for Encode)
  let path = d3.path();
  const r = read_radius;
  let gap = 0;
  let index = 0;
  while (read.is_unit[index] == 0 && index < read.is_unit.length) {
    gap = read.contig[index];
    index += 1;
  }
  // Current ID of the contig
  let contig = read.contig[index];
  let current_unit = read.unit[index];
  let start = start_pos[contig] - Math.PI / 2;
  let radian = start + bp_scale(unit_length * current_unit);
  //if (gap != 0) {
  //   path.moveTo(
  //     contig_radius * Math.cos(radian),
  //     contig_radius * Math.sin(radian)
  //   );
  //   path.lineTo(read_radius * Math.cos(radian), read_radius * Math.sin(radian));
  // } else {
  path.moveTo(read_radius * Math.cos(radian), read_radius * Math.sin(radian));
  //  }
  gap = 0;
  for (; index < read.unit.length; index += 1) {
    if (read.is_unit[index] == 0) {
      gap = read.contig[index];
      continue;
    }
    const diff = Math.abs(read.unit[index] - current_unit);
    const is_end_to_start =
      (current_unit < 30 && read.unit[index] > ends[read.contig[index]] - 30) ||
      (current_unit > ends[read.contig[index]] - 30 && read.unit[index] < 30);
    current_unit = read.unit[index];
    if ((read.contig[index] == contig && diff < 100) || is_end_to_start) {
      radian = start + bp_scale(unit_length * read.unit[index]);
      path.lineTo(r * Math.cos(radian), r * Math.sin(radian));
    } else {
      // Change contig. Connect them.
      const new_radian = start_pos[read.contig[index]];
      radian =
        new_radian + bp_scale(unit_length * read.unit[index]) - Math.PI / 2;
      // Bezier Curve to new point from here.
      const control_radius =
        handle_points[contig][read.contig[index]] - Math.PI / 2;
      const control_x = handle_points_radius * Math.cos(control_radius);
      const control_y = handle_points_radius * Math.sin(control_radius);
      contig = read.contig[index];
      start = start_pos[contig] - Math.PI / 2;
      path.quadraticCurveTo(
        control_x,
        control_y,
        r * Math.cos(radian),
        r * Math.sin(radian)
      );
    }
  }
  // if (gap != 0) {
  //   path.moveTo(
  //     contig_radius * Math.cos(radian),
  //     contig_radius * Math.sin(radian)
  //   );
  //   path.lineTo(read_radius * Math.cos(radian), read_radius * Math.sin(radian));
  // }
  return path.toString();
};

const getNumOfGapRead = (reads) => {
  // Input: [JSON object]
  // Output: Num
  // Requirements: each element should be 'read' object.
  // Return numbers of reads which is just Gap.
  return reads.every((read) => read.is_unit == 0).length;
};

const criticalpairToPath = (
  cp,
  handle_points,
  bp_scale,
  start_pos,
  unit_length
) => {
  const r = read_radius;
  let path = d3.path();
  // Move to contig1
  const contig1 = cp["contig1"];
  const contig1_start_angle = start_pos[contig1["contig"]] - Math.PI / 2;
  const start_angle_1 =
    contig1_start_angle + bp_scale(unit_length * contig1["start_unit"]);
  const end_angle_1 =
    contig1_start_angle + bp_scale(unit_length * contig1["end_unit"]);
  path.moveTo(r * Math.cos(start_angle_1), r * Math.sin(start_angle_1));
  path.arc(0, 0, r, start_angle_1, end_angle_1);
  // Bezier Curve to contig2.
  const contig2 = cp["contig2"];
  const contig2_start_angle = start_pos[contig2["contig"]] - Math.PI / 2;
  const start_angle_2 =
    contig2_start_angle + bp_scale(unit_length * contig2["start_unit"]);
  const control_radius =
    handle_points[contig1["contig"]][contig2["contig"]] - Math.PI / 2;
  const control_x = handle_points_radius * Math.cos(control_radius);
  const control_y = handle_points_radius * Math.sin(control_radius);
  path.quadraticCurveTo(
    control_x,
    control_y,
    r * Math.cos(start_angle_2),
    r * Math.sin(start_angle_2)
  );
  const end_angle_2 =
    contig2_start_angle + bp_scale(unit_length * contig2["end_unit"]);
  path.arc(0, 0, r, start_angle_2, end_angle_2);
  path.quadraticCurveTo(
    control_x,
    control_y,
    r * Math.cos(start_angle_1),
    r * Math.sin(start_angle_1)
  );
  return path.toString();
};

const confluentregionToPath = (
  cr,
  handle_points,
  bp_scale,
  start_pos,
  unit_length
) => {
  const to_r = read_radius + 50;
  const over_r = read_radius + 60;
  let path = d3.path();
  const contig = cr["pos"];
  const contig_start_angle = start_pos[contig["contig"]] - Math.PI / 2;
  const start_angle =
    contig_start_angle +
    bp_scale(unit_length * contig["start_unit"]) -
    confluent_margin;
  const end_angle =
    contig_start_angle +
    bp_scale(unit_length * contig["end_unit"]) +
    confluent_margin;
  if (contig["direction"] == "UpStream") {
    const overshoot = start_angle - Math.PI / 70;
    path.moveTo(
      read_radius * Math.cos(start_angle),
      read_radius * Math.sin(start_angle)
    );
    path.lineTo(to_r * Math.cos(start_angle), to_r * Math.sin(start_angle));
    path.arc(0, 0, to_r, start_angle, overshoot, true);
    path.lineTo(over_r * Math.cos(overshoot), over_r * Math.sin(overshoot));
    path.arc(0, 0, over_r, overshoot, end_angle, false);
    path.lineTo(
      read_radius * Math.cos(end_angle),
      read_radius * Math.sin(end_angle)
    );
    path.closePath();
  } else if (contig["direction"] == "DownStream") {
    const overshoot = end_angle + Math.PI / 70;
    path.moveTo(
      read_radius * Math.cos(start_angle),
      read_radius * Math.sin(start_angle)
    );
    path.lineTo(over_r * Math.cos(start_angle), over_r * Math.sin(start_angle));
    path.arc(0, 0, over_r, start_angle, overshoot, false);
    path.lineTo(to_r * Math.cos(overshoot), to_r * Math.sin(overshoot));
    path.arc(0, 0, to_r, overshoot, end_angle, true);
    path.lineTo(
      read_radius * Math.cos(end_angle),
      read_radius * Math.sin(end_angle)
    );
    path.closePath();
  }
  return path.toString();
};

const crToPath = (cr, handle_points, bp_scale, start_pos, unit_length) => {
  // Input: JSON object, JSON object, Integer
  // Output: String
  // Requirements: Critical region object, scales
  // Return the path btw critical region, or confluent path.
  if (cr.hasOwnProperty("CP")) {
    return criticalpairToPath(
      cr["CP"],
      handle_points,
      bp_scale,
      start_pos,
      unit_length
    );
  } else if (cr.hasOwnProperty("CR")) {
    return confluentregionToPath(
      cr["CR"],
      handle_points,
      bp_scale,
      start_pos,
      unit_length
    );
  } else {
    console.log(`Error ${cr}`);
    return 1;
  }
};

const htgap = (read) => {
  let sum = 0;
  if (read.is_unit[read.is_unit.length - 1] == 0) {
    sum += read.contig[read.is_unit.length - 1];
  }
  if (read.is_unit[0] == 0) {
    sum += read.contig[0];
  }
  return sum;
};

const calcGap = (reads) => {
  const len = reads.length;
  const sum = reads.map((read) => htgap(read)).reduce((acc, x) => acc + x, 0);
  return (sum / len) * 2;
};

const kFormatter = (num) => {
  return Math.abs(num) > 999
    ? Math.sign(num) * (Math.abs(num) / 1000).toFixed(1) + "k"
    : Math.sign(num) * Math.abs(num);
};

const contigToHTML = (contig) => {
  const start = kFormatter(contig["start_unit"] * 100);
  const end = kFormatter(contig["end_unit"] * 100);
  const direction = contig["direction"];
  return `<div><ul>
<li>Start:${start} bp</li>
<li>End:${end} bp</li>
<li>Direction:${direction} </li>
</ul>
</div>`;
};

const criticalpairToHTML = (cp, idx) => {
  const header = `<div class = critical-region><div>CP:${idx}</div>`;
  const contig1 = contigToHTML(cp["contig1"]);
  const contig2 = contigToHTML(cp["contig2"]);
  return header + contig1 + contig2 + "</div>";
};

const confluentregionToHTML = (cr, idx) => {
  const header = `<div class = critical-region><div>CR:${idx}</div>`;
  const contig = contigToHTML(cr["pos"]);
  return header + contig + "</div>";
};

const crToHTML = (cr, cluster) => {
  // Input: JSON object, Array
  // Output: String
  // Requirements: Critical region object
  // Return the HTML contents corresponds to the given cr.
  if (cr.hasOwnProperty("CP")) {
    return criticalpairToHTML(cr["CP"], cluster);
  } else if (cr.hasOwnProperty("CR")) {
    return confluentregionToHTML(cr["CR"], cluster);
  } else {
    console.log(`Error ${cr}`);
    return "Error";
  }
};

const calcCoverageOf = (reads, contigs, members) => {
  // Input: List of JSON object, List of JSON object, JSON object.
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
  let results = contigs.map((covs) => {
    const len = covs.length;
    let coverage = Array.from({ length: len }, (_) => 0);
    return { id: covs.id, length: covs.length, cov: coverage };
  });
  for (const read of reads) {
    for (let idx = 0; idx < read.is_unit.length; idx += 1) {
      if (read.is_unit[idx] == 1) {
        const c = read.contig[idx];
        const p = read.unit[idx];
        results[c].cov[p] += 1;
      }
    }
  }
  const positions = members.flatMap((member) => {
    const cr = member.cr;
    if (cr.hasOwnProperty("CP")) {
      return [cr.CP.contig1, cr.CP.contig2];
    } else if (cr.hasOwnProperty("CR")) {
      return [cr.CR.pos];
    }
  });
  const MAX_COV = 10000;
  for (const pos of positions) {
    const contig = pos.contig;
    const start = pos.start_unit;
    results[contig].cov[start] = MAX_COV;
  }
  return results.map((data) => {
    data.cov = data.cov
      .map((d, idx) => {
        return { cov: d, pos: idx };
      })
      .filter((d) => d.cov > 1)
      .map((d) => {
        if (d.cov === MAX_COV) {
          d.cov = 0;
        }
        return d;
      });
    return data;
  });
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

const plotData = (dataset, repeats, unit_length) =>
  Promise.all([dataset, repeats].map((file) => d3.json(file)))
    .then(([values, repeats]) => {
      // Unpack
      // This is array.
      const contigs = values.contigs;
      // This is also an array.
      // const reads = values.reads;
      // Or select reads as you like.
      const reads = values.reads.map(convertRead);
      const clusters = values.clusters;
      // Calculate coordinate.
      const bp_scale = calcScale(contigs);
      const coverage_scale = calcCovScale(contigs);
      const start_pos = calcStartPosition(contigs);
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
      // Draw contigs.
      const selection_on_each_contig = contigs_layer
        .selectAll(".contig")
        .data(contigs)
        .enter()
        .append("g")
        .attr("transform", (contig) => `rotate(${start_pos[contig.id]})`);
      selection_on_each_contig
        .append("path")
        .attr("class", "contig")
        .attr("id", (c) => `contig-${c.id}`)
        .attr("d", (contig) => {
          const end = bp_scale(contig.length);
          const arc = d3
            .arc()
            .innerRadius(contig_radius - contig_thick)
            .outerRadius(contig_radius)
            .startAngle(0)
            .endAngle(end);
          return arc();
        })
        .attr("fill", (c) => d3.schemeCategory10[c.id % 10])
        .attr("opacity", 0.2);
      selection_on_each_contig
        .append("text")
        .attr("dy", -30)
        .attr("dx", 20)
        .append("textPath")
        .attr("href", (c) => `#contig-${c.id}`)
        .text((c) => `${c.name}`);
      const selection_on_ticks = selection_on_each_contig
        .append("g")
        .selectAll(".tick")
        .data((contig) => {
          const len = contig.length / reference_tick;
          return Array.from({ length: len }).map((_, idx) => {
            return {
              angle: bp_scale(idx * reference_tick),
              position: idx * reference_tick,
            };
          });
        })
        .enter()
        .append("g")
        .attr("class", "tick")
        .attr(
          "transform",
          (d) =>
            `rotate(${
              (180 * d.angle) / Math.PI + 180
            }) translate(0,${contig_radius})`
        );
      selection_on_ticks.append("line").attr("y2", 10).attr("stroke", "black");
      selection_on_ticks
        .filter((d) => d.position % reference_break === 0 && d.position != 0)
        .append("text")
        .attr("transform", (d) => {
          if (d.angle > Math.PI) {
            return `rotate(-90) translate(-60,0)`;
          } else {
            return `rotate(90) translate(20,0)`;
          }
        })
        .text((d) => `${kFormatter(d.position)}`);
      // Draw repeat.
      contigs_layer
        .selectAll(".repeats")
        .data(repeats.flatMap((rp) => rp.reps))
        .enter()
        .append("path")
        .attr("class", "repeats")
        .attr("d", (repeat) => {
          const arc = d3
            .arc()
            .innerRadius(contig_radius - contig_thick - 3)
            .outerRadius(contig_radius + 3)
            .startAngle(start_pos[repeat.id] + bp_scale(repeat.start))
            .endAngle(start_pos[repeat.id] + bp_scale(repeat.end));
          return arc();
        })
        .attr("fill", "gray");
      // Draw coverage
      const selection_on_each_coverage = coverage_layer
        .selectAll(".coverage")
        .data(contigs)
        .enter()
        .append("g")
        .attr("class", (contig) => `coverage-${contig.id}`)
        .attr(
          "transfrom",
          (contig) => `rotate(${(start_pos[contig.id] * 180) / Math.PI})`
        );
      selection_on_each_coverage
        .append("path")
        .attr("class", "coverage")
        .attr("d", (contig) => {
          const arc = d3
            .lineRadial()
            .angle((_, i) => bp_scale(i * unit_length))
            .radius((d) => coverage_scale(d));
          return arc(contig.coverages);
        })
        .attr("fill", "none")
        .attr("stroke", (c) => d3.schemeCategory10[c.id % 10]);
      selection_on_each_coverage.each((contig) => {
        const max_cov = Math.max(...contig.coverages);
        const domain = [0, max_cov];
        const range = [-coverage_min, -coverage_scale(max_cov)];
        const local_scale = d3.scaleLinear().domain(domain).range(range);
        d3.selectAll(`.coverage-${contig.id}`).call(
          d3.axisLeft(local_scale).tickFormat(d3.format(".2s")).ticks(2)
        );
      });
      // Draw start/stop reads.
      const selection_on_each_start_stop = start_stop_layer
        .selectAll(".start-stop")
        .data(contigs)
        .enter()
        .append("g")
        .attr("class", (contig) => `start-stop-${contig.id}`)
        .attr(
          "transform",
          (contig) => `rotate(${(start_pos[contig.id] * 180) / Math.PI})`
        );
      selection_on_each_start_stop
        .append("path")
        .attr("class", "start-stop")
        .attr("d", (contig) => {
          const arc = d3
            .lineRadial()
            .angle((_, i) => bp_scale(i * unit_length))
            .radius((d) => readnum_scale(d));
          return arc(contig.start_stop);
        })
        .attr("fill", "none")
        .attr("stroke", (c) => d3.schemeCategory10[c.id % 10]);
      selection_on_each_start_stop.each((contig) => {
        const max_num = Math.max(...contig.start_stop);
        const domain = [0, max_num];
        const range = [-start_stop_radius, -readnum_scale(max_num)];
        const local_scale = d3.scaleLinear().domain(domain).range(range);
        const tV = Array.from({
          length: Math.floor(max_num / start_stop_break),
        }).map((_, idx) => (idx + 1) * start_stop_break);
        d3.selectAll(`.start-stop-${contig.id}`).call(
          d3.axisLeft(local_scale).tickValues(tV)
        );
      });
      const max_lengths = contigs.map((c) => c.length / unit_length);
      // Draw reads
      read_layer
        .selectAll(".read")
        .data(reads)
        .enter()
        .append("path")
        .attr("class", "read")
        .attr("d", (read) =>
          readToPath(
            read,
            handle_points,
            bp_scale,
            start_pos,
            unit_length,
            max_lengths
          )
        )
        .attr("fill", "none")
        .attr("opacity", 0.1)
        .attr("stroke", (read) => "black");
      // Draw critical regions.
      const max_cluster_id = Math.max(...clusters.map((d) => d.id));
      let is_active = Array.from({ length: max_cluster_id }, (_) => false);
      const critical_regions = clusters.flatMap((d) => d.members);
      cr_layer
        .selectAll(".cr")
        .data(critical_regions)
        .enter()
        .append("path")
        .attr("class", "cr")
        .attr("d", (d) =>
          crToPath(d.cr, handle_points, bp_scale, start_pos, unit_length)
        )
        .attr("stroke", (d) => d3.schemeCategory10[(d.cluster + 1) % 10])
        .attr("stroke-width", (member) =>
          member.cr.hasOwnProperty("CP") ? 5 : 5
        )
        .attr("opacity", (member) =>
          member.cr.hasOwnProperty("CP") ? 0.4 : 0.5
        )
        .attr(
          "fill",
          (member) => d3.schemeCategory10[(member.cluster + 1) % 10]
        )
        .on("click", function (member) {
          // Check if the cluster is already clicked.
          const cluster = member.cluster;
          const active = is_active[cluster];
          if (active) {
            is_active[cluster] = false;
            temp_coverage_layer.selectAll(`.cluster-${cluster}`).remove();
            cr_info.select(`.cluster-${cluster}`).remove();
          } else {
            is_active[cluster] = true;
            const supporting_reads = reads.filter((r) => r.cluster == cluster);
            const members = clusters.find((cl) => cl.id === cluster).members;
            const coverages = calcCoverageOf(
              supporting_reads,
              contigs,
              members
            );
            const [argmax, max] = coverages
              .map((cov) =>
                cov.cov.reduce((y, x) => (x.cov > y.cov ? x : y), {
                  cov: 0,
                  pos: 0,
                })
              )
              .reduce(
                ([contig, y], x, c) => (x.cov > y.cov ? [c, x] : [contig, y]),
                [0, { cov: 0, pos: 0 }]
              );
            const total_coverage = contigs.find((c) => c.id == argmax)
              .coverages[max.pos];
            console.log(argmax, max, total_coverage);
            let contents = critical_regions.reduce((acc, member, idx) => {
              if (member.cluster == cluster) {
                const html = crToHTML(member.cr, idx);
                acc += html;
              }
              return acc;
            }, `<div>Cluster ${cluster}</div>`);
            const count = supporting_reads.length;
            const meangap = calcGap(supporting_reads);
            const support = `<div>Supporing Reads:${count}</div>`;
            const cov_ratio =
              `<div>Fraction:${max.cov}/${total_coverage} ` +
              `(${kFormatter(max.pos * unit_length)} bp)</div>`;
            const gap = `<div>Mean gap length:${meangap.toFixed(1)}</div>`;
            contents += support + gap + cov_ratio;
            const info_tip = cr_info
              .insert("div", ":first-child")
              .classed(`cluster-${cluster} cluster-parent`, true);
            info_tip.append("div").html(contents);
            info_tip
              .append("div")
              .attr("class", "info-tip-clustercolor")
              .append("svg")
              .attr("width", 200)
              .attr("height", 20)
              .append("rect")
              .attr("x", 0)
              .attr("y", 0)
              .attr("width", 200)
              .attr("height", 20)
              .attr("rx", 2)
              .attr("fill", d3.schemeCategory10[(cluster + 1) % 10]);
            console.log(members);
            temp_coverage_layer
              .selectAll(".tempcoverage")
              .data(coverages)
              .enter()
              .append("path")
              .attr("class", `cluster-${cluster}`)
              .attr("d", (coverage) => {
                const start = start_pos[coverage.id];
                const arc = d3
                  .lineRadial()
                  .angle((d) => start + bp_scale(d.pos * unit_length))
                  .radius((d) => coverage_scale(d.cov))
                  .defined((d, i, data) => {
                    if (i + 1 === data.length) {
                      return true;
                    } else {
                      return data[i].pos + 1 === data[i + 1].pos;
                    }
                  });
                return arc(coverage.cov);
              })
              .attr("fill", "none")
              .attr("stroke-width", 2)
              .attr("stroke", d3.schemeCategory10[(cluster + 1) % 10]);
          }
        });
      info
        .append("div")
        .attr("class", "numofgapread")
        .append("p")
        .text(`Gap Read:${getNumOfGapRead(reads)} out of ${reads.length}`);
      return scales;
    })
    .then(
      (ok) => ok,
      (why) => console.log(why)
    );
