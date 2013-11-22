var bluePinkColorscale = d3.scale.linear()
         .domain([-3, 0, 3])
         .range(["blue", "white", "red"]);


NetNodeModel = Backbone.Model.extend({
   urlRoot : "/network/nodes",
   idAttribute : "netfig_node_id"
});

NetNodeView = Backbone.View.extend({
    tagname : "div",
    className : "node",
    render : function() {
        this.$el.html("<a class='nodelabel' target=_search href='/search?q=" + this.model.get("node_obj_id") + "'>" + this.model.get("node_obj_id") + "</a>");
        this.$el.css("left", this.model.get("node_pos_left"));
        this.$el.css("top", this.model.get("node_pos_top"));
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

// model for basic figure info without actual data
NetFigureMetadata = Backbone.Model.extend({

});

NetFigureMetadataList = Backbone.Collection.extend({
    model: NetFigureMetadata,
    url: "/network/figures"
});

NetFigureMetadataView = Backbone.View.extend({
    tagName: "option",
    render: function() {
        this.el.innerHTML = this.model.get("netfig_id").toString() + ":  " + this.model.get("title");
        $(this.el).attr("value", this.model.get("netfig_id") );
        return this;
    }
});

NetFigureMetadataListView = Backbone.View.extend({
    tagName: "select",
    render: function() {
      var that = this;
      $(this.el).empty()
      this.collection.each( function(f) {
          var nfv = new NetFigureMetadataView({ model: f });
          $(that.el).append(nfv.render().el);
        });
    }
})

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
        var that = this;
        $("#figtitle").text(this.model.get("title"));
        $("#figdescription").text(this.model.get("description"));

        var endpointOptions = {
                    isSource:true,
                    isTarget: true,
                    connectorStyle : { strokeStyle:"#666", curviness: 100 },
                    paintStyle: { radius: 4, fillStyle: 'black'}
        };

        var newNodeView = function(node) {
            var nv = new NetNodeView({ model : node, parentGraphView : this });
	        node_element = nv.render().el;
            //alert(node);
            that.$el.children("#" + node.get("node_compartment")).append( node_element );

            jsPlumb.draggable(node_element, { containment:"parent",
                stop :
                    function(event, ui) {
//                        var d = nv.model.get("data");
//                        d["position"]
                        nv.model.set("node_pos_left", ui.position.left);
                        nv.model.set("node_pos_top", ui.position.top);
                        nv.model.save();

                    }
            } );

            var endpoint = jsPlumb.addEndpoint(node_element, endpointOptions);
            return nv;
        }

        that.child_nodeviews = {};
        this.model.nodes.each(function(node) {
            that.child_nodeviews[node.get("netfig_node_id")] = newNodeView(node);
        } , this);

        that.child_edgeviews = {};
        this.model.edges.each( function(edge) {
           //  var ev = new NetEdgeView({model: edge, parentGraphView : this});
           // ev.render();

           var overlayType;
           var paintStyle;
           if (edge.get("edge_type") === "activates") {
               overlayType = [ [ "Arrow", { location:0.9, fillStyle: "green" } ] ];
               paintStyle = { strokeStyle:"darkgreen", lineWidth:1 };
           }
           else if (edge.get("edge_type") === "inhibits") {
               overlayType = [ [ "Diamond", { location:0.9, fillStyle: "red", width: 10, length: 10} ] ];
               paintStyle = { strokeStyle:"darkred", lineWidth: 1 };
           }
           else {
               overlayType = [];
               paintStyle = { strokeStyle:"#666", lineWidth:1 };
           }

           jsPlumb.connect( {
               source: that.child_nodeviews[edge.get("source_node_id")].el,
               target : that.child_nodeviews[edge.get("target_node_id")].el,
               paintStyle: paintStyle,
               endpointStyle:{ radius: 4, fillStyle: '#666' },
               anchors:[ [ "Perimeter", { shape:"Rectangle" } ], ["Perimeter", { shape:"Rectangle" } ]],
               overlays: overlayType
               }
           );

        }, this);

    },


    setColors: function(contrastSet) {
        // load expression data
        var that = this;

        $.getJSON("/network/expression/" + this.model.get("netfig_id").toString(),
            function(data)  {
                // iterate through nodes, setting color if it matches
                _.each(that.child_nodeviews, function(node_view) {
                    if (node_view.model.get("node_obj_idtype") === "entrez_gene_symbol") {
                    var gene_symbol = node_view.model.get("node_obj_id");
                    $(node_view.el).css("background", bluePinkColorscale( data[gene_symbol][contrastSet] ));
                    }
            });
        }
        );

    }
});


// var $neighbourList = $("<div class='selectionWindow'></div>").text(" pick here.. ");
// $("body").append($neighbourList);
// $neighbourList.draggable().resizable();
