// This is a tiny library to visualize multipartite structure of mitochondrial genome of plant.
// There is only one function exposed as a public API, namely `plotData`.
// It takes three arguments and autoamtically plot the data to a DOM(default:<div id="plot"></div>).  The size of the canvas is 1000x1000 and some other DOM elements would be overwritten(<div id="info"></div>).

// The the API is two URIs to serialized JSON objects, `dataset` and `repeats`, and one integer value named `unit_length`.

// The `dataset` value should consist of three JSON objects named `contigs`, `reads`, and `clusters`. Below I explain how these objects should be structured:

// - `contigs`: contigs is an array of contig object. A contig object is a JSON-object having elements as follows:
// ```json
// {
// name:String,
// id: Integer,
// length: Integer,
// coverages:Array<Integer>,
// start_stop<Integer>,
// }
// ```
// - `reads`: reads is an array of read. Each read is a JSON-object as follwos:
// ```json
// {
// name:String,
// units:Array<Unit>,
// cluster:Integer,
// }
// ```
// and Unit is either
// ```json
// {
// G:Integer,
// }
// ```
// or
// ```json
// {
// E:[Interger;2]
// }
// ```

// - `clusters`: clusters is an array of array of CriticalRegion. Each CriticalRegion is either
// ```
// {
// "CP" -> ContigPair,
// }
// ```
// or
// ```
// {
// "CR" -> ConfluentRegion,
// }
// ```

// - ContigPair:
// ```
// {
// contig1:Position,
// contig2:Position,
// reads:HashSet<String>,
// }
// ```
// - ConfluentRegion
// ```
// {
// pos:Position,
// reads:HashSet<String>,
// }
// ```

// - Position
// ```
// {
// contig:Integer,
// start_unit:Integer,
// end_unit:Integer,
// length:Integer,
// direction: Either UpStream or DownStream.
// }
// ```

const plotData = (dataset, alignments, unit_length) =>
  Promise.all([dataset, alignments].map((file) => d3.json(file)))
    .then(([values, contig_alns]) => {
      const width = 1000;
      const height = 400 + values.contigs.length * 100;
      // Margin between contigs
      const left_margin = 100;
      const left_contig_margin = 100 + left_margin;
      const contig_margin = 100;
      const contig_thick = 10;
      const coverage_thick = 5;
      const read_thick = 4;
      const eplison = 5.001;
      const confluent_margin = 0.01;
      // The maximum contig length.
      const max_contig_length = 600;
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
      const gap_scale = d3
        .scaleLog()
        .domain([gap_min_size, gap_max_size])
        .range([gap_min, gap_max])
        .clamp(true);
      const tick_top_margin = 20;
      const tick_margin = 50;
      const master_circle = 20;
      const master_circle_margin = 25;
      const mc_total_margin = master_circle + master_circle_margin;
      const colorScheme = d3.schemePaired;

      const svg = d3
        .select("#plot")
        .append("svg")
        .attr("width", width)
        .attr("height", height);
      const contigs_layer = svg
        .append("g")
        .attr("transform", `translate(${left_contig_margin},${height / 2})`)
        .attr("class", "contigs");
      const coverage_layer = svg
        .append("g")
        .attr("transform", `translate(${left_contig_margin},${height / 2})`)
        .attr("class", "coverages");
      const start_stop_layer = svg
        .append("g")
        .attr("transform", `translate(${left_contig_margin},${height / 2})`)
        .attr("class", "start-stop-read");
      const read_layer = svg
        .append("g")
        .attr("transform", `translate(${left_contig_margin},${height / 2})`)
        .attr("class", "read");
      const aln_layer = svg
        .append("g")
        .attr("transform", `translate(${left_contig_margin},${height / 2})`)
        .attr("class", "critical-region");

      const tooltip = d3
        .select("body")
        .append("div")
        .attr("class", "tooltip")
        .style("opacity", 0);
      const info = d3.select("#info");

      const calcScale = (contigs) => {
        // Input: Array of JSON object
        // Output: d3 Scale object
        // Requirements: each input object should have "length" attribute
        // Convert base pair into radian
        const max = Math.max(...contigs.map((c) => c.length));
        return d3.scaleLinear().domain([0, max]).range([0, max_contig_length]);
      };

      const calcStartPosition = (contigs) => {
        // Input: Array of JSON object.
        // Output: Array[Num]
        // Requirements: each input object should have "length" attribute
        // Map from the index(id) of contig into the start position of the contig(the y-axis).
        // Note that the start position may be negative, as the center is (left_margin, height/2)
        const center = ((contigs.length - 1) * contig_margin) / 2;
        let start_pos = Array.from({ length: contigs.length }, (_) => 0);
        contigs.forEach((contig, idx) => {
          start_pos[contig.id] = idx * contig_margin - center;
        });
        return start_pos;
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
        return d3
          .scaleLinear()
          .domain([0, max])
          .range([coverage_min, coverage_max]);
      };

      const calcReadNumScale = (contigs) => {
        // Input: Array on JSON object
        // Output: d3.scale object
        // Requirements: Each object in the argument should have an array of integer, which is
        // named "start_stop."
        // Calculate the scale for start/stop vizualization.
        const total = contigs
          .flatMap((c) => c.start_stop)
          .reduce((x, y) => x + y);
        const num = contigs
          .map((c) => c.start_stop.length)
          .reduce((x, y) => x + y);
        const max = Math.max(...contigs.flatMap((c) => c.start_stop));
        return d3
          .scaleLog()
          .domain([min_read_num, max])
          .range([min_radius, max_radius])
          .clamp(true);
      };

      const readToPath = (
        read,
        handle_points,
        bp_scale,
        start_pos,
        unit_length
      ) => {
        // Input: JSON object, Array[Array[Num]], d3.scale, Array[Num], Num
        // Output: String
        // Requirements: read should have units attribute, each of which elements
        // should have either "G"(for Gap) or "E"(for Encode)
        let path = d3.path();
        let units = Array.from(read.units).reverse();
        let gap = 0;
        let unit = {};
        while (!unit.hasOwnProperty("E")) {
          unit = units.pop();
          if (unit == undefined) {
            return "";
          } else if (unit.hasOwnProperty("G")) {
            gap = unit.G;
          }
        }
        // Current ID of the contig
        let contig = unit.E[0];
        let current_unit = unit.E[1];
        let y = start_pos[contig];
        let x = bp_scale(unit_length * unit.E[1]);
        if (gap != 0) {
          path.moveTo(x, y + gap_scale(gap));
          path.lineTo(x, y);
        } else {
          path.moveTo(x, y);
        }
        for (unit of units.reverse()) {
          if (unit.hasOwnProperty("G")) {
            continue;
          }
          const diff = Math.abs(unit.E[1] - current_unit);
          if (unit.E[0] == contig && diff < 50) {
            // x =  bp_scale(unit_length*unit.E[1]);
            // path.lineTo(x, y);
          } else {
            x = bp_scale(unit_length * current_unit);
            path.lineTo(x, y);
            // Change contig. Connect them.
            const new_y = start_pos[unit.E[0]];
            const new_x = bp_scale(unit_length * unit.E[1]);
            // Bezier Curve to new point from here.
            const control_x = (x + new_x) / 2;
            const control_y = handle_points[contig][unit.E[0]];
            contig = unit.E[0];
            path.quadraticCurveTo(control_x, control_y, new_x, new_y);
            x = new_x;
            y = new_y;
          }
          current_unit = unit.E[1];
        }
        return path.toString();
      };

      const selectRead = (read, unitlen) => {
        // Input: JSON object, Num
        // Output: boolean
        // Requirements: input object should have "units" property,
        // which is actually vector of object with "Gap" or "Encode" property.
        // Filter read as you like.
        return true;
      };

      const getNumOfGapRead = (reads) => {
        // Input: [JSON object]
        // Output: Num
        // Requirements: each element should be 'read' object.
        // Return numbers of reads which is just Gap.
        return reads.filter((read) => {
          let units = Array.from(read.units);
          let unit = {};
          while (!unit.hasOwnProperty("E")) {
            unit = units.pop();
            if (unit == undefined) {
              return true;
            }
          }
          return false;
        }).length;
      };

      const htgap = (read) => {
        let sum = 0;
        if (read.units[read.units.length - 1].hasOwnProperty("G")) {
          sum += read.units[read.units.length - 1].G;
        }
        if (read.units[0].hasOwnProperty("G")) {
          sum += read.units[0].G;
        }
        return sum;
      };

      const calcGap = (reads) => {
        const len = reads.length;
        const sum = reads
          .map((read) => htgap(read))
          .reduce((acc, x) => acc + x, 0);
        return (sum / len) * 2;
      };

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

      // Unpack
      // This is array.
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
      const contigs = values.contigs.sort(compare_name);
      for (const c of contigs) {
        console.log(c.name, c.id, c.length);
      }
      // This is also an array.
      // const reads = values.reads;
      // Or select reads as you like.
      const reads = values.reads.filter((r) => selectRead(r, unit_length));
      // const critical_regions = [values.critical_regions[selected_region]];
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
      // Draw coverage
      coverage_layer
        .selectAll(".coverage")
        .data(contigs)
        .enter()
        .append("path")
        .attr("class", "coverage")
        .attr("d", (contig) => {
          const height = start_pos[contig.id];
          const path = d3
            .line()
            .x((_, i) => bp_scale(i * unit_length))
            .y((cov) => height - coverage_scale(cov));
          return path(contig.coverages);
        })
        .attr("fill", "none")
        .attr("stroke", (c) => d3.schemeCategory10[c.id % 10]);
      // Draw reads
      read_layer
        .selectAll(".read")
        .data(reads)
        .enter()
        .append("path")
        .attr("class", "read")
        .attr("d", (read) =>
          readToPath(read, handle_points, bp_scale, start_pos, unit_length)
        )
        .attr("fill", "none")
        .attr("opacity", 0.3)
        .attr("stroke", (read) => "black");
      info
        .append("div")
        .attr("class", "numofgapread")
        .append("p")
        .text(`Gap Read:${getNumOfGapRead(reads)} out of ${reads.length}`);
      // Draw contigs
      const each_contig_group = contigs_layer
        .selectAll(".contig")
        .data(contigs)
        .enter()
        .append("g")
        .attr("transform", (c) => `translate(0, ${start_pos[c.id]})`);
      each_contig_group
        .append("path")
        .attr("class", "contig")
        .attr("d", (contig) => {
          const end = bp_scale(contig.length);
          let path = d3.path();
          path.moveTo(0, 0);
          path.lineTo(end, 0);
          path.lineTo(end, -contig_thick);
          path.lineTo(0, -contig_thick);
          path.closePath();
          return path.toString();
        })
        .attr("fill", "gray")
        .attr("opacity", 0.2);
      each_contig_group
        .append("text")
        .attr("class", "contig")
        .attr("text-anchor", "end")
        .text((contig) => contig.name);
      // Draw ticks.
      const b_tick = svg
        .append("g")
        .attr("class", "scale")
        .attr("transform", `translate(0,${tick_top_margin})`);
      b_tick.append("text").text("Base Pair Scale");
      {
        b_tick
          .append("g")
          .attr("transform", "translate(50,5)")
          .call(d3.axisBottom(bp_scale).tickFormat(d3.format(".2s")).ticks(4));
      }
      const p_tick = svg
        .append("g")
        .attr("class", "scale")
        .attr("transform", `translate(0,${tick_margin + tick_top_margin})`);
      p_tick.append("text").text("Position of the master circle.");
      {
        const gradient = p_tick
          .selectAll("contigGradient")
          .data(contig_alns.references)
          .enter()
          .append("defs")
          .attr("class", "contigGradient")
          .append("linearGradient")
          .attr("id", (_, idx) => `contig-${idx}`)
          .attr("x1", "0%")
          .attr("x2", "100%")
          .attr("y1", "0%")
          .attr("y2", "0%");
        gradient
          .append("stop")
          .attr("offset", "0%")
          .attr("stop-color", (_, idx) => colorScheme[(2 * idx + 1) % 10]);
        gradient
          .append("stop")
          .attr("offset", "100%")
          .attr("stop-color", (_, idx) => colorScheme[(2 * idx + 2) % 10]);
        const p_group = p_tick
          .selectAll("masterCircle")
          .data(contig_alns.references)
          .enter()
          .append("g")
          .attr(
            "transform",
            (_, i) =>
              `translate(50,${master_circle + i * master_circle_margin})`
          )
          .attr("class", "masterCircle");
        p_group
          .append("text")
          .attr("class", "masterCircleName")
          .text((c) => c.id);
        p_group
          .append("rect")
          .attr("transform", `translate(0,5)`)
          .attr("x", 0)
          .attr("width", (c) => bp_scale(c.length))
          .attr("y", 0)
          .attr("height", (c) => 10)
          .attr("fill", (_, idx) => `url(#contig-${idx})`);
      }
      const occupied = contig_alns.references.length * mc_total_margin;
      const c_tick = svg
        .append("g")
        .attr("class", "scale")
        .attr(
          "transform",
          `translate(0,${tick_margin + 2 * tick_top_margin + occupied})`
        );
      c_tick.append("text").text("Coverage Scale");
      {
        const cscale = d3
          .scaleLinear()
          .domain([0, 500])
          .range([0, coverage_scale(500) - coverage_scale(0)]);
        c_tick
          .append("g")
          .attr("transform", `translate(50,5)`)
          .call(d3.axisBottom(cscale).tickFormat(d3.format(".2s")).ticks(3));
      }
      const g_tick = svg
        .append("g")
        .attr("class", "scale")
        .attr(
          "transform",
          `translate(0,${2 * tick_margin + 3 * tick_top_margin + occupied})`
        );
      g_tick.append("text").text("Gap Scale");
      {
        const gscale = d3
          .scaleLog()
          .domain([gap_min, 5 * gap_max])
          .range([0, 5 * (gap_max - gap_min)]);
        g_tick
          .append("g")
          .attr("transform", `translate(50,5)`)
          .call(d3.axisBottom(gscale).tickFormat(d3.format(".2s")).ticks(1));
      }
      // Draw alignments
      const gradient = aln_layer
        .selectAll("alnGradient")
        .data(contig_alns.alignments)
        .enter()
        .append("defs")
        .attr("class", "alnGadient")
        .append("linearGradient")
        .attr("id", (_, idx) => `align-${idx}`)
        .attr("x1", (aln) => {
          const rfr = contig_alns.references.find(
            (r) => r.id == aln.reference_name
          );
          if (aln.query_start < aln.query_end) {
            const fraction = aln.ref_start / (aln.ref_end - aln.ref_start);
            return `${-100 * fraction}%`;
          } else {
            const fraction = (rfr.length - aln.ref_end) / rfr.length;
            return `${100 * (1 + fraction)}%`;
          }
        })
        .attr("x2", (aln) => {
          const rfr = contig_alns.references.find(
            (r) => r.id == aln.reference_name
          );
          if (aln.query_start < aln.query_end) {
            const fraction = (rfr.length - aln.ref_end) / rfr.length;
            return `${100 * (1 + fraction)}%`;
          } else {
            const fraction = aln.ref_start / (aln.ref_end - aln.ref_start);
            return `${-100 * fraction}%`;
          }
        })
        .attr("y1", "0%")
        .attr("y2", "0%");
      gradient
        .append("stop")
        .attr("offset", "0%")
        .attr("stop-color", (aln) => {
          const color = contig_alns.references
            .map((refr, idx) => [refr.id, idx])
            .find((d) => d[0] == aln.reference_name);
          return colorScheme[(2 * color[1] + 1) % 10];
        });
      gradient
        .append("stop")
        .attr("offset", "100%")
        .attr("stop-color", (aln) => {
          const color = contig_alns.references
            .map((refr, idx) => [refr.id, idx])
            .find((d) => d[0] == aln.reference_name);
          return colorScheme[(2 * color[1] + 2) % 10];
        });
      aln_layer
        .selectAll(".alignments")
        .data(contig_alns.alignments)
        .enter()
        .append("rect")
        .attr("x", (aln) => {
          const start = bp_scale(Math.min(aln.query_start, aln.query_end));
          return start;
        })
        .attr("y", (aln) => {
          const elm = contigs.find((c) => c.name === aln.query_name);
          return start_pos[elm.id] - contig_thick;
        })
        .attr("height", contig_thick)
        .attr("width", (aln) => {
          if (aln.query_start < aln.query_end) {
            return bp_scale(aln.query_end - aln.query_start);
          } else {
            return bp_scale(aln.query_start - aln.query_end);
          }
        })
        .attr("fill", (_, idx) => `url(#align-${idx})`);
      return scales;
    })
    .then(
      (ok) => ok,
      (why) => console.log(why)
    );
