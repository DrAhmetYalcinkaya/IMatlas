const instance = tippy("#graph", {
	followCursor: true,
	arrow: false,
	theme: "blank",
	hideOnClick: false
});

var graph = document.getElementById("graph");
graph.addEventListener("wheel", function(event) {
	Shiny.setInputValue("zooming", {"deltaY" : event.deltaY,
		"x" : event.x,
		"y" : event.y}, 
	{priority: "event"});
	return true;
}, true);

Shiny.addCustomMessageHandler("tooltip", function(name) {
	instance[0].setContent(name);
});

Shiny.addCustomMessageHandler("debug", function(text) {
    console.log(text);
});

$(document).on('shiny:sessioninitialized', function(event) {
    var height = window.innerHeight;
    console.log(height);
    Shiny.setInputValue("screenHeight", height);
});

Shiny.addCustomMessageHandler("cursor", function(cursor){
    document.getElementById("graph").style.cursor = cursor;
});




