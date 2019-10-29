const width = 1000;
const height = 1000;
const MAX_POWER=0.05;
const MIN_POWER=0.02;
const RADIUS = 200;
const svg = d3.select("#plot")
      .append("svg")
      .attr("width",width)
      .attr("height",height);
const graph_layer = svg
      .append("g")
      .attr("transform",`translate(${height/2},${width/2})`);

const info = d3.select("#info");

const remove_dedupe = (values)=>{
    let dedupe = new Set();
    for (const d of values){
        dedupe.add(d.pos1);
        dedupe.add(d.pos2);
    }
    return Array.from(dedupe).map(p => {
        return {p:p};
    });
};

const drag = simulation => {
    function dragstarted(d) {
        if (!d3.event.active) simulation.alphaTarget(0.3).restart();
        d.fx = d.x;
        d.fy = d.y;
    }
    function dragged(d) {
        d.fx = d3.event.x;
        d.fy = d3.event.y;
    }
    function dragended(d) {
        if (!d3.event.active) simulation.alphaTarget(0);
        d.fx = null;
        d.fy = null;
    }
    return d3.drag()
        .on("start", dragstarted)
        .on("drag", dragged)
        .on("end", dragended);
};

const plotData = (dataset) =>
      Promise.all([dataset]
                  .map(file => d3.tsv(file)))
      .then(([vs]) => {
          const values = vs.filter(d => d.pos1 != d.pos2)
                .filter(d => d.mlpvalue > 20);
          console.log(values.slice(1,10));
          console.log(values.length);
          const scale = d3.interpolateCividis;
          const points = remove_dedupe(values);
          console.log(points.length);
          console.log(points.slice(0,10));
          const links = values.map(d => {
              return {
                  source: d.pos1,
                  target: d.pos2,
                  strength: d.mlpvalue
              };
          });
          const max_pos = Math.max(...points.map(d => d.p));
          const radian_scale = d3.scaleLinear()
                .domain([0, max_pos])
                .range([0, 2*Math.PI]);
          const nodes = points.map(d => {
              const rad = radian_scale(d.p);
              return {
                  p:d.p,
                  x: Math.cos(rad) * RADIUS,
                  y: Math.sin(rad) * RADIUS,
              };
          });
          const node_layer = graph_layer
                .append("g")
                .attr("class","nodes")
                .selectAll(".node")
                .data(nodes)
                .enter()
                .append("circle")
                .attr("class","node")
                .attr("fill",d => scale(d.p/max_pos))
                .attr("r",3);
          // const link_layer = graph_layer
          //       .append("g")
          //       .attr("class","links")
          //       .selectAll(".link")
          //       .data(links)
          //       .enter()
          //       .append("line")
          //       .attr("class","link");
          const max = Math.max(...links.map(l => l.strength));
          const power_scale = d3.scaleLinear()
                .domain([0,max])
                .range([MIN_POWER,MAX_POWER]);
          const simulation = d3.forceSimulation()
                .nodes(nodes)
                .force("charge",d3.forceManyBody()
                       .strength(() => -0.03))
                .force("link",
                       d3.forceLink()
                       .links(links)
                       .id(d => d.p)
                       .strength(l => power_scale(l.strength))
                      )
                .on("tick", () => {
                    // link_layer
                    //     .attr("x1", l => l.source.x)
                    //     .attr("y1", l => l.source.y)
                    //     .attr("x2", l => l.target.x)
                    //     .attr("y2", l => l.target.y);
                    node_layer.attr("cx", d => d.x)
                        .attr("cy", d => d.y);
                });
          node_layer.call(drag(simulation));
      })
      .then(ok => ok,
            why => console.log(why));

