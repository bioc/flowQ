
var bg="#7896b1";
var fg="#c9d8e5";

function toggleImage(id){
    var hide = "none";
    var show = "inline";
    var obj = document.getElementById("img_"+id);
    var state = obj.style.display;
    if(state == show){
	obj.style.display=hide;
    }else{
	obj.style.display=show;
    }
}


function toggleDetails(iRange, jRange, k, tobj, tid){
    var hide = "none";
    var show = "block";
    for(var i=1; i<=iRange; i++){
	for(var j=1; j<=jRange; j++){
	    var obj = document.getElementById("button_"+k+"_"+i+"_"+j);
	    var hobj = document.getElementById(tid);
	    var row1 = document.getElementById("row_"+k+"_"+i+"_"+j+"_1");
	    var row2 = document.getElementById("row_"+k+"_"+i+"_"+j+"_2");
	    var image = document.getElementById("img_"+k+"_"+i+"_"+j);
	    var state = obj.style.display;
	    if(state == show){
		obj.style.display=hide;
		image.style.display=hide;
		tobj.style.display=hide;
		hobj.style.display=show;
		row1.style.width="0px";
		row2.style.width="0px";
	    }else{
		obj.style.display=show;
		tobj.style.display=hide;
		hobj.style.display=show;
		row1.style.width="0px";
		row2.style.width="0px";
	    }
	}
    }
}



