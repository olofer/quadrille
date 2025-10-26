const DonothingDemo = {

    NAME: "(none)",

    init: function () {
        console.log("Initializing: " + this.NAME);
    },

    reset: function () {
    },

    evolve: function (tsim_, dt_) {
    },

    draw: function (ctx_, tsim_, xmin_, xmax_, ymin_, ymax_) {
    },

    print_stats: function (ctx_, tsim_) {
        //ctx_.fillText("demo: " + this.NAME, 5, 35);
    },

    handle_key_down: function (e) {
        const code = e.keyCode;
        const key = e.key;
    },
}
