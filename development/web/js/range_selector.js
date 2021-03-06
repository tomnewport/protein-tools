// Range Selector plugin (slider) by Tom Newport
// Corrected using JSHint
// Call using slider( element ) passing element as a selector, HTML element or JS element.


function slider( slider_element ){
		var slider_container = $(slider_element);
		slider_container.addClass("range_selector")

		var slider_color = slider_container.data("color");

		var slider_span = 	$("	<div class='span rs_container'></div>");
		var slider_pre  = $("<span class='pre'></span>");
		var slider_mid  = $("<span class='mid'></span>");
		var slider_post = $("<span class='post'></span>");

		slider_pre.css('background-color', slider_color);
		slider_mid.css('background-color', slider_color);
		slider_post.css('background-color', slider_color);

		var slider_start = 	$("<div class='slider start'><div class='inner'></div></div>");
		var slider_end = 	$("<div class='slider end'><div class='inner'></div></div>");

		slider_start.css('background-color', slider_color);
		slider_end.css('background-color', slider_color);

		var label_start = 	$("<input class='rs_label start' value ='10'/>");
		var label_end = 	$("<input class='rs_label end' value = '20'/>");

		label_start.css('color', slider_color);
		label_end.css('color', slider_color);

		slider_span.append(slider_pre, slider_mid, slider_post, slider_start, slider_end, label_start, label_end);

		slider_container.append(slider_span);

		slider_container.data('start_handle', slider_start );
		slider_container.data('end_handle'  , slider_end );
		slider_container.data('start_label' , label_start );
		slider_container.data('end_label'   , label_end );
		slider_container.data('pre_span'    , slider_pre  );
		slider_container.data('mid_span'    , slider_mid  );
		slider_container.data('post_span'   , slider_post );

		slider_start.data(	'parent', slider_container);
		slider_end.data(	'parent', slider_container);
		label_start.data(	'parent', slider_container);
		label_end.data(		'parent', slider_container);
		slider_pre.data(	'parent', slider_container);
		slider_mid.data(	'parent', slider_container);
		slider_post.data(	'parent', slider_container);

				label_start.change( function(e){
					var tgt = $(e.target).data('parent');
					var v = parseInt($(e.target).val());
					tgt.data('start', v);
					tgt.attr('data-start', v);
					slider_update( tgt );
				});

				label_end.change( function(e){
					var tgt = $(e.target).data('parent');
					var v = parseInt($(e.target).val());
					tgt.data('end', v);
					tgt.attr('data-end', v);
					slider_update( tgt );
				});

				//0 = no drag, 1 = start, 2 = end
				slider_container.data('dragmode', 0);

				//Change drag mode:
				slider_start.mousedown(function(e){
					var active = $(e.currentTarget).data('parent');
					active.data('dragmode', 1);
				});
				slider_start.mouseup(function(e){
					var active = $(e.currentTarget).data('parent');
					active.data('dragmode', 0);
				});
				slider_end.mousedown(function(e){
					var active = $(e.currentTarget).data('parent');
					active.data('dragmode', 2);
				});
				slider_end.mouseup(function(e){
					var active = $(e.currentTarget).data('parent');
					active.data('dragmode', 0);
				});

				slider_mid.mousedown(function(e){
					var active = $(e.currentTarget).data('parent');
					active.data('dragmode', 3);
					active.data('dragstart', -1);
					$(e.currentTarget).addClass("moving");
				});

				slider_mid.mouseup(function(e){
					var active = $(e.currentTarget).data('parent');
					active.data('dragmode', 0);
					$(e.currentTarget).removeClass("moving");
				});

				slider_container.mouseleave(function(e){
					var active = $(e.currentTarget);
					active.data('dragmode', 0);
				});

				slider_container.mousemove(function(e){
					var tgt = $(e.currentTarget);
					if (tgt.data('dragmode') !== 0){
						
						var cx = e.clientX;
						var tx = tgt.find('.span.rs_container').offset().left;
						var tw = tgt.find('.span.rs_container').innerWidth();
						var rp = Math.max(0, Math.min((cx - tx)/tw, 1));

						var smin = tgt.data('min');
						var smax = tgt.data('max');
						var aap = smin + (rp * (smax - smin));

						if (tgt.data('int') === true){
							aap = Math.round(aap);
						}
						if (tgt.data('dragmode') == 1){
							tgt.data('start', aap);
							tgt.attr('data-start', aap);
						} else if (tgt.data('dragmode') == 2){
							tgt.data('end', aap);
							tgt.attr('data-end', aap);
						} else if (tgt.data('dragmode') == 3){
							if (tgt.data('dragstart') == -1){
								tgt.data('dragstart', aap);
								tgt.data('startstart', tgt.data('start'));
								tgt.data('endstart', tgt.data('end'));
							}
							var dragdiff = aap - tgt.data('dragstart');

							var st = Math.max(tgt.data('startstart') + dragdiff, tgt.data('min'));
							var en = Math.min(tgt.data('endstart') + dragdiff, tgt.data('max'));

							tgt.data('start', st);
							tgt.attr('data-start', st);

							tgt.data('end', en);
							tgt.attr('data-end', en);
							
							//tgt.data('end', aap);
							//tgt.attr('data-end', aap);
						} 
						slider_update( tgt );
					}
				});
				slider_update( slider_container );
				return slider_container;
		}
		function slider_update(slider){
			//Update Slider
			var pfunc = parseFloat;

			if (slider.data('int') === true){
				pfunc = parseInt;
			}
				var smax =   pfunc(slider.attr('data-max'));
				var smin =   pfunc(slider.attr('data-min'));
				var sstart = pfunc(slider.attr('data-start'));
				var send =   pfunc(slider.attr('data-end'));

			//Enforce logic:

			//sstart cannot be greater than send

			sstart = Math.min(send - 1, sstart);

			//sstart cannot be less than smin

			sstart = Math.max(sstart, smin);

			//send cannot be less than sstart

			send = Math.max(sstart + 1, send);

			//send cannot be greater than smax

			send = Math.min(smax, send);

			//Return checked logic to variables
			slider.attr('data-max', smax);
			slider.attr('data-min', smin);
			slider.attr('data-start', sstart);
			slider.attr('data-end', send);
			slider.data('max', smax);
			slider.data('min', smin);
			slider.data('start', sstart);
			slider.data('end', send);

			slider.data('start_label').val(sstart);
			slider.data('end_label').val(send);
			//Move labels
			var start_perc = 100 * ((sstart - smin) / (smax - smin));
			var end_perc = 100 * ((send - smin) / (smax - smin));
			slider.data('start_handle').css('right', 'calc('+(100-start_perc)+'% - 4px)');
			slider.data('end_handle').css('left', 'calc('+end_perc+'% - 4px)');
			slider.data('start_label').css('right', (100-start_perc)+'%');
			slider.data('end_label').css('left', end_perc+'%');

			slider.data('pre_span').css('width', start_perc + '%');
			slider.data('mid_span').css('width', (end_perc - start_perc) + '%');
			slider.data('post_span').css('width', (100 - end_perc) + '%');

			slider.trigger( "update" );

		}