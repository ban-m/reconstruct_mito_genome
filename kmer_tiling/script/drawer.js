const width = 1000;
const height = 1000;
const top_margin = 50;
const left_margin = 50;
const scale_width = 500;
const scale_height = 50;
const radius = 20;
const font_size = 10;
const size = 7;
const svg = d3.select(".plot")
      .append("svg")
      .attr("width",width)
      .attr("height",height);
const colorScheme = d3.interpolateViridis;
const canvas = svg.append("g").attr("transform",`translate(${width/2},${height/2})`);
const color_bar = svg.append("g")
      .attr("transform",`translate(${top_margin},${left_margin})`);

const maxOfList = list => 
      list.reduce((acc,current) => Math.max(acc,current),0);
const sumOfList = list => list.reduce((acc,current) => acc + current,0);

const countToArray = count => {
    const dist = count.dist;
    const len = count.count.length;
    return count.count.map((c,i) => {
        return {"count":c,
                "dist":dist,
                "rad":Math.PI * i * 2 / len,
               };});
};


const plotKmer = (kmer,ave_count, x, y) => {
    canvas
        .append("g")
        .attr("class","kmer")
        .selectAll(".observe")
        .data(kmer.dists.flatMap(countToArray))
        .enter()
        .append("circle")
        .attr("class","observe")
        .attr("cx", d => x + radius * d.dist * Math.cos(d.rad))
        .attr("cy", d => y + radius * d.dist * Math.sin(d.rad))
        .attr("r", size)
        .attr("fill", d => colorScheme(Math.min(d.count/ave_count,1)));
};


const plotColorScheme = ave_count => {
        const len = 100;
    const colorSample = Array.from({"length":len},(v,i) => i)
          .map( i => colorScheme(i * ave_count / len ));
    const scale = d3.scaleLinear()
          .domain([0,len])
          .range([0,scale_width]);
    const bar_width = scale_width/len - 1;// 1 px margin
    color_bar.selectAll("colorBar")
        .data(colorSample)
        .enter()
        .append("rect")
        .attr("x", (c,idx) => scale(idx))
        .attr("y", 0)
        .attr("width", bar_width)
        .attr("height", scale_height)
        .attr("fill", c => c);
    color_bar
        .append("text")
        .attr("x",0)
        .attr("y",scale_height + font_size)
        .text("0 count");
    color_bar
        .append("text")
        .attr("x",scale_width - 3 * font_size)
        .attr("y",scale_height + font_size)
        .text(`${ave_count.toFixed(2)} count`);

};

const plotJsonFile = data => {
    const total_length = data.reduce((acc,kmer) => acc + kmer.dists.reduce((acc, dist) => acc + dist.count.length,0),0);
    const ave_count = sumOfList(data.map((kmer) => sumOfList(kmer.dists.map( dist => sumOfList(dist.count))))) / total_length;
    console.log(`Ave count:${ave_count}`);
    for (const kmer of data) {
        const x = (Math.random()-1/2) * width;
        const y = (Math.random()-1/2) * height;
        plotKmer(kmer, ave_count,x,y);
    }
    plotColorScheme(ave_count);
};

const plotConsectiveKmer = data => {
    console.log(data);
    const total_length = data.reduce((acc,kmer) => acc + kmer.dists.reduce((acc, dist) => acc + dist.count.length,0),0);
    const ave_count = sumOfList(data.map((kmer) => sumOfList(kmer.dists.map( dist => sumOfList(dist.count))))) / total_length;
    console.log(`Ave count:${ave_count}`);
    const margin = Math.max(height,width)/6;
    const len = data.length;
    const direction = Math.PI /4;
    const xinc = radius * Math.sin(direction);
    const yinc = -radius * Math.cos(direction);
    const xstart = width/6;
    const ystart = height/6;
    const trajectory = Array.from({"length":len},(v,i) => i)
          .map(i => {
              return {"x":xstart + xinc * i,
                      "y":ystart + yinc * i};
          });
    data.forEach((kmer,idx) =>  {
        const {x,y} = trajectory[idx];
        plotKmer(kmer, ave_count,x,y);
    });
    plotColorScheme(ave_count);

};

const plotConsectiveData = (jsonfile) =>
      d3.json(jsonfile)
      .then(plotConsectiveKmer)
      .then(ok => console.log("OK"),
            ng => console.log(ng));

const plotScatterData = (jsonfile) =>
      d3.json(jsonfile)
      .then(plotJsonFile)
      .then(ok => console.log("OK"),
            ng => console.log(ng));
