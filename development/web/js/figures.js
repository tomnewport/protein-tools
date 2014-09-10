//Javascript library to plot figures from JSON files.
//By Tom Newport 2014

function figures_overview(element, data){
	fig_oview = $(element);
	fig_oview.addClass("fig_overview");
	fig_svg = $("<div><svg width='100%' height='100%' viewBox='0 0 1600 150'></svg></div>");
	fig_oview.append(fig_svg);
	s = Snap( fig_svg.find("svg")[0] );

	annotations = data["annotations"];
	l("Drawing " + annotations.length + " annotations.");
	x = 0;
	for (aid in annotations){
		a = annotations[aid];
		span = s.path("M"+ x + " 75L" + (1600*(a.start/aa_length)) + " 75");
		span.attr({"stroke": "#fff", "stroke-opacity": 1, "stroke-width": 4, "vector-effect": "non-scaling-stroke"});
		a_svg = s.rect(1600*(a.start/aa_length), 10, 1600*((a.end - a.start)/aa_length), 130, 10,10  );
		a_svg.data("annotation",a);
		a_svg.attr({"cursor":"pointer","fill": "#ffffff", "fill-opacity": 0.5, "stroke": "#fff", "stroke-opacity": 1, "stroke-width": 4, "vector-effect": "non-scaling-stroke" });
		a_svg.hover(
		function(e){
			Snap(e.target).animate({"fill-opacity": 1}, 200);
			d = Snap(e.target).data("annotation");
		}, 
		function(e){
			Snap(e.target).animate({"fill-opacity": 0.5}, 200);
		});
		a_svg.click(
		function(e){
			alert(Snap(e.target).data("annotation").name);
		});
		x = 1600*(a.end/aa_length);
	}
	span = s.path("M"+ x + " 75L" + 1600 + " 75");
	span.attr({"stroke": "#fff", "stroke-opacity": 1, "stroke-width": 4, "vector-effect": "non-scaling-stroke"});
	return fig_oview;
}