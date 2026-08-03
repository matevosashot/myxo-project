import glob
import os
import time

import h5py
import plotly.graph_objects as go


data_dir = "/home/ashmat/Projects/myxo-project/data/scale_change_4.1/"
slider_key = "lM"
filename = "scale_change_P"


t0 = time.perf_counter()

hdf5_files = glob.glob(os.path.join(data_dir, "*.h5"))
data = []
for hdf5_file in hdf5_files:
    with h5py.File(hdf5_file, "r") as f:
        datum = {}
        metadata = dict(f.attrs.items())
        datum["metadata"] = metadata

        datum["P_contrib_diag"] = f["P_contrib_diag"][:]
        datum["Q_contrib_diag"] = f["Q_contrib_diag"][:]
        datum["grid_1d"] = f["grid_1d"][:]
        data.append(datum)


P_max = max(d["P_contrib_diag"].max() for d in data)
P_min = min(d["P_contrib_diag"].min() for d in data)
Q_max = max(d["Q_contrib_diag"].max() for d in data)
Q_min = min(d["Q_contrib_diag"].min() for d in data)

t1 = time.perf_counter()
print(f"[timing] load data: {t1 - t0:.3f}s")

##########################
### edit below


# sort data so the slider sweeps through a meaningful parameter
sorted_data = sorted(data, key=lambda d: d["metadata"][slider_key])
keys = [d["metadata"][slider_key] for d in sorted_data]
n = len(sorted_data)

# build traces: P (left scene) first, then Q (right scene)
#
# Convention used here: z is indexed as z[ix, iy], i.e. the 0-th axis of z
# corresponds to x and the 1-st axis corresponds to y.
# Plotly's go.Surface places z[i, j] at (x=x[j], y=y[i]), so we transpose
# z to match the (ix, iy) convention.
left_traces = [
    go.Surface(
        x=datum["grid_1d"],
        y=datum["grid_1d"],
        z=datum["P_contrib_diag"].T,
        colorscale="viridis",
        visible=(i == 0),
        showscale=False,
        scene="scene",
        name=f"P {slider_key}={keys[i]:.3g}",
    )
    for i, datum in enumerate(sorted_data)
]
right_traces = [
    go.Surface(
        x=datum["grid_1d"],
        y=datum["grid_1d"],
        z=datum["Q_contrib_diag"].T,
        colorscale="viridis",
        visible=(i == 0),
        showscale=False,
        scene="scene2",
        name=f"Q {slider_key}={keys[i]:.3g}",
    )
    for i, datum in enumerate(sorted_data)
]
traces = left_traces + right_traces
t2 = time.perf_counter()
print(f"[timing] build {len(traces)} surface traces: {t2 - t1:.3f}s")


def title_for(i):
    meta = sorted_data[i]["metadata"]
    return (
        f"{slider_key} = {meta[slider_key]:.3g}  |  "
        f"lm = {meta['lm']:.3g}  |  "
        f"lM = {meta['lM']:.3g}"
    )


# one slider step per dataset, toggling visibility for both P and Q traces.
steps = []
for i, datum in enumerate(sorted_data):
    visibility = (
        [j == i for j in range(n)]   # left (P)
        + [j == i for j in range(n)] # right (Q)
    )
    steps.append(
        dict(
            method="update",
            args=[
                {"visible": visibility},
                {"title.text": title_for(i)},
            ],
            label=f"{keys[i]:.3g}",
        )
    )

sliders = [
    dict(
        active=0,
        currentvalue={
            "prefix": f"{slider_key} = ",
            "visible": True,
            "font": {"size": 16, "color": "black"},
            "xanchor": "left",
        },
        pad={"t": 50},
        steps=steps,
    )
]
t3 = time.perf_counter()
print(f"[timing] build {len(steps)} slider steps: {t3 - t2:.3f}s")


# build figure with two scenes side-by-side
fig = go.Figure(data=traces)
common_axes = dict(
    xaxis_title="x",
    yaxis_title="y",
    # xaxis=dict(range=[-40, 40]),
    # yaxis=dict(range=[-40, 40]),
)
fig.update_layout(
    sliders=sliders,
    scene=dict(
        domain={"x": [0.0, 0.48], "y": [0.0, 1.0]},
        zaxis_title="P contribution",
        # zaxis=dict(range=[P_min, P_max]),
        **common_axes,
    ),
    scene2=dict(
        domain={"x": [0.52, 1.0], "y": [0.0, 1.0]},
        zaxis_title="Q contribution",
        # zaxis=dict(range=[Q_min, Q_max]),
        **common_axes,
    ),
    title=title_for(0),
    annotations=[
        dict(
            text="P_contrib_diag",
            x=0.24, y=1.0, xref="paper", yref="paper",
            showarrow=False, font=dict(size=14),
        ),
        dict(
            text="Q_contrib_diag",
            x=0.76, y=1.0, xref="paper", yref="paper",
            showarrow=False, font=dict(size=14),
        ),
    ],
)
t4 = time.perf_counter()
print(f"[timing] create figure + update_layout: {t4 - t3:.3f}s")


# Lock cameras between the two 3D scenes: dragging/zooming/panning one
# mirrors the camera onto the other.
sync_camera_js = r"""
var gd = document.getElementById('{plot_id}');
if (gd) {
    var syncing = false;
    gd.on('plotly_relayout', function(ed) {
        if (syncing || !ed) return;
        var update = null;
        if (ed['scene.camera']) {
            update = {'scene2.camera': ed['scene.camera']};
        } else if (ed['scene2.camera']) {
            update = {'scene.camera': ed['scene2.camera']};
        }
        if (update) {
            syncing = true;
            Plotly.relayout(gd, update).then(function() { syncing = false; })
                .catch(function() { syncing = false; });
        }
    });
}
"""

script_dir = os.path.dirname(os.path.abspath(__file__))
fig.write_html(
    os.path.join(script_dir, filename + ".html"),
    post_script=sync_camera_js,
)
t5 = time.perf_counter()
print(f"[timing] write_html: {t5 - t4:.3f}s")

print(f"[timing] total: {t5 - t0:.3f}s")
