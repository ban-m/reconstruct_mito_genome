pub const TEMPLATE: &str = r#"<!DOCTYPE html>
<meta charset="utf-8">
<html>
  <head>
    <script src="https://d3js.org/d3.v5.min.js"></script>
    <link rel="stylesheet" type="text/css" href="/viewer/style.css">
  </head>
  <body>
    <div class = "title">
      Circos plot of Plant mitochondria.
    </div>
    <div id = "plot"></div>
    <div id = "info"></div>
    <script src="/script/circos.js"></script>
    <script>
      const dataset = "/viewer/data.json";
      const repeats = "/viewer/repeats.json";
      const unit_length = 100;
      plotData(dataset,repeats,unit_length);
    </script>
  </body>
</html>
"#;

pub const TEMPLATE_LINEAR: &str = r#"<!DOCTYPE html>
<meta charset="utf-8">
<html>
  <head>
    <script src="https://d3js.org/d3.v5.min.js"></script>
    <link rel="stylesheet" type="text/css" href="/viewer/style.css">
  </head>
  <body>
    <div class = "title">
      Circos plot of Plant mitochondria.
    </div>
    <div id = "plot"></div>
    <div id = "info"></div>
    <script src="/script/linear.js"></script>
    <script>
      const dataset = "/viewer/read_data.json";
      const repeats = "/viewer/contig_alns.json";
      const unit_length = 100;
      plotData(dataset,repeats,unit_length);
    </script>
  </body>
</html>
"#;


pub const STYLE: &str = r#".title{
    font-family: sans-serif;
    font-size: large;
}

.scale{
    font-family: sans-serif;
}

.scale text{
    font-size: medium;
}

.tick{
    stroke-width:1;
}

.contig{
    opacity:0.5;
}

.repeats{
    opacity:0.4;
}


.numofgapread{
    font-family: sans-serif;
    font-size: medium;
}

.critical_region{
}

.tooltip {
    position: absolute;			
    text-align: left;
    padding: 2px;				
    font: 14px sans-serif;		
    background: lightsteelblue;	
    border: 0px;		
    border-radius: 8px;			
    pointer-events: none;			
}
"#;
