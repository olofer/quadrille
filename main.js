/*

Demo of the planeview.js grid plotter with an interactive animation.

Hot-swap which demo is executing by pressing the (shift) tab key.

Each "demo" object must have these member functions:
  init, reset, draw, evolve, print_stats, handle_key_down

*/

const LIST_OF_DEMOS = [BouyancyDemo, GravityDemo, DonothingDemo];

for (let i = 0; i < LIST_OF_DEMOS.length; i++) {
    LIST_OF_DEMOS[i].init();
}

let demoIndex = 0;

const DEFAULT_BBOX = [-5.0, 5.0, -3.0, 3.0];

const canvas = document.getElementById("canvas");
let PV = createDefaultPlaneView(canvas.width, canvas.height);

let viewZoom = 1.00;
let viewEta = 0.125;
let bbox = [DEFAULT_BBOX[0], DEFAULT_BBOX[1], DEFAULT_BBOX[2], DEFAULT_BBOX[3]];

const FPS = 50.0;  // 20ms refresh intervals
const dt = 0.0020;  // 2.0ms simulator time-stepping

let filtered_fpsval = 0.0;
const filter_beta = 0.990;

let tsim = 0.0;
let frame = 0;
let startTime = Date.now();

function resizeWindow() {
    canvas.width = window.innerWidth - 24;
    canvas.height = window.innerHeight - 24;
    PV = createDefaultPlaneView(canvas.width, canvas.height);
}

function reset_state() {
    tsim = 0.0;
    frame = 0;
    startTime = Date.now();
    LIST_OF_DEMOS[demoIndex].reset();
}

function printMasterStats(ctx, timestamp, fpsval) {
    ctx.font = "17px Arial";
    ctx.fillStyle = "white";
    ctx.fillText("demo: " + LIST_OF_DEMOS[demoIndex].NAME +
        " -- time: " + timestamp.toFixed(3) + " sec" +
        " [" + fpsval.toFixed(1) + " fps]", 5, 15);
}

function refresh() {

    let currentTime = Date.now();
    let elapsedTime = currentTime - startTime;
    startTime = currentTime;

    if (elapsedTime <= 0.0) return;

    let elapsedSec = elapsedTime / 1000.0;

    while (elapsedSec > 0.0) {

        LIST_OF_DEMOS[demoIndex].evolve(tsim, dt);

        tsim += dt;
        elapsedSec -= dt;
    }

    PV.autoZoom(bbox[0], bbox[1], bbox[2], bbox[3], viewZoom, viewEta);

    const ctx = canvas.getContext("2d");

    PV.drawGrid(ctx);
    PV.setTransform(ctx);

    LIST_OF_DEMOS[demoIndex].draw(ctx, tsim, PV.xmin, PV.xmax, PV.ymin, PV.ymax);

    PV.unitTransform(ctx);
    filtered_fpsval = filter_beta * filtered_fpsval + (1.0 - filter_beta) * (1000.0 / elapsedTime);
    printMasterStats(ctx, tsim, filtered_fpsval);
    LIST_OF_DEMOS[demoIndex].print_stats(ctx, tsim);
    PV.setTransform(ctx);

    frame++;
}

function keyDownEvent(e) {
    let code = e.keyCode;
    let key = e.key;

    if (key == "Tab") {
        if (e.shiftKey) demoIndex--; else demoIndex++;
        if (demoIndex == LIST_OF_DEMOS.length) demoIndex = 0;
        if (demoIndex == -1) demoIndex = LIST_OF_DEMOS.length - 1;
        // reset_state();
        e.preventDefault();
        return;
    }

    if (key == 'r' || key == 'R') {
        reset_state();
        return;
    }

    if (key == 'z' || key == 'Z') { // reset view
        bbox = [DEFAULT_BBOX[0], DEFAULT_BBOX[1], DEFAULT_BBOX[2], DEFAULT_BBOX[3]];
        viewZoom = 1.0;
        return;
    }

    if (key == 'c' || key == 'C') { // center view (do not reset zoom level)
        bbox = [DEFAULT_BBOX[0], DEFAULT_BBOX[1], DEFAULT_BBOX[2], DEFAULT_BBOX[3]];
        return;
    }

    if (code == 38 && e.shiftKey) {
        viewZoom *= 0.80;
        return;
    }

    if (code == 40 && e.shiftKey) {
        viewZoom *= 1.25;
        return;
    }

    const dx = (bbox[1] - bbox[0]) / 12.0;
    const dy = (bbox[3] - bbox[2]) / 12.0;
    const dmove = viewZoom * (dx + dy) / 2;

    if (code === 39) // right
    {
        bbox[0] += dmove;
        bbox[1] += dmove;
        return;
    }

    if (code === 37) // left
    {
        bbox[0] -= dmove;
        bbox[1] -= dmove;
        return;
    }

    if (code === 38) // up
    {
        bbox[2] += dmove;
        bbox[3] += dmove;
        return;
    }

    if (code === 40) // down
    {
        bbox[2] -= dmove;
        bbox[3] -= dmove;
        return;
    }

    LIST_OF_DEMOS[demoIndex].handle_key_down(e);
}

function keyUpEvent(e) {
    let code = e.keyCode;
    let key = e.key;
}

resizeWindow();
reset_state();

window.addEventListener('resize', resizeWindow);
window.addEventListener('keydown', keyDownEvent);
window.addEventListener('keyup', keyUpEvent);

let refresher_id = setInterval(refresh, 1000 / FPS);
