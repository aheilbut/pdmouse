
NetNodeModel = Backbone.Model.extend({
});

NetNodeView = Backbone.View.extend({
    tagname : "div",
    className : "node",
    render : function() {
        this.$el.html(this.model.get("node_obj_id"));
        return this;
    }
});


NetEdgeModel = Backbone.Model.extend({
});

NetNodeList = Backbone.Collection.extend({
    model : NetNodeModel
});

NetEdgeList = Backbone.Collection.extend({
    model: NetEdgeModel
});


 NetFigure = Backbone.Model.extend({
    urlRoot : "/network/figure",

    parse: function(response, options) {
        this.nodes = new NetNodeList(response.nodes);
        this.edges = new NetEdgeList(response.edges);
        return response;

//        this.nodes.reset( this.get("nodes"));
 //       this.edges.reset( this.get("edges"));
    },

    addNode: function(nodeType, nodeId) {

    },

    addEdge: function(nodeType, nodeId ) {

    }
});

 NetFigureView = Backbone.View.extend({
    el : "#netfig",

    render : function() {
        $(this.el).find("#figtitle").text(this.model.get("title"));
        $(this.el).find("#figdescription").text(this.model.get("description"));

        this.model.nodes.each( function(node) {
            var nv = new NetNodeView({ model : node, parentGraphView : this });
	         node_element = nv.render().el;
            //alert(node);
            this.$el.children("#" + node.get("node_compartment")).append( node_element );

            jsPlumb.draggable(node_element, { containment:"parent" });

            /*	    jsPlumb.addEndpoint(node_element, { endpoint : "Rectangle",
                                    isSource: true,
                                    maxConnections: -1,
                                    isTarget: true,
                                    hoverPaintStyle:{ fillStyle:"#449999" },
                                    hoverClass: "endpointHoverClass",
                                    paintStyle: { fillStyle: "black",
                                              width : 15,
                                              height: 15 } } );
            */
        }, this);


    }
});


// var $neighbourList = $("<div class='selectionWindow'></div>").text(" pick here.. ");
// $("body").append($neighbourList);
// $neighbourList.draggable().resizable();
