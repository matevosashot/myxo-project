import bisect
import glob
import os
import time

import h5py
import plotly.graph_objects as go


slider_key = "L"
filename = "space_change_normalized_P"


t0 = time.perf_counter()

data_dir = "/home/ashmat/Projects/myxo-project/data/space_change/"
hdf5_files = glob.glob(os.path.join(data_dir, "*.h5"))
data = []
for hdf5_file in hdf5_files:
    with h5py.File(hdf5_file, "r") as f:
        datum = {}
        metadata = dict(f.attrs.items())
        datum["metadata"] = metadata
        datum["P_contrib_diag"] = f["P_contrib_diag"][:]
        datum["grid_1d"] = f["grid_1d"][:]
        data.append(datum)


data_dir = "/home/ashmat/Projects/myxo-project/data/space_change_normalized/"
hdf5_files = glob.glob(os.path.join(data_dir, "*.h5"))
data_normalized = []
for hdf5_file in hdf5_files:
    with h5py.File(hdf5_file, "r") as f:
        datum = {}
        metadata = dict(f.attrs.items())
        datum["metadata"] = metadata
        datum["P_contrib_diag"] = f["P_contrib_diag"][:]
        datum["grid_1d"] = f["grid_1d"][:]
        data_normalized.append(datum)


z_max_left = max(d["P_contrib_diag"].max() for d in data)
z_min_left = min(d["P_contrib_diag"].min() for d in data)
z_max_right = max(d["P_contrib_diag"].max() for d in data_normalized)
z_min_right = min(d["P_contrib_diag"].min() for d in data_normalized)

t1 = time.perf_counter()
print(f"[timing] load data: {t1 - t0:.3f}s")

##########################
### edit below


# sort each dataset by the slider key
sorted_left = sorted(data, key=lambda d: d["metadata"][slider_key])
sorted_right = sorted(data_normalized, key=lambda d: d["metadata"][slider_key])
keys_left = [d["metadata"][slider_key] for d in sorted_left]
keys_right = [d["metadata"][slider_key] for d in sorted_right]
n_left = len(sorted_left)
n_right = len(sorted_right)

# build traces: left scene first, then right scene
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
        name=f"raw {slider_key}={keys_left[i]:.3g}",
    )
    for i, datum in enumerate(sorted_left)
]
right_traces = [
    go.Surface(
        x=datum["grid_1d"],
        y=datum["grid_1d"],
        z=datum["P_contrib_diag"].T,
        colorscale="viridis",
        visible=(i == 0),
        showscale=False,
        scene="scene2",
        name=f"norm {slider_key}={keys_right[i]:.3g}",
    )
    for i, datum in enumerate(sorted_right)
]
traces = left_traces + right_traces
t2 = time.perf_counter()
print(f"[timing] build {len(traces)} surface traces: {t2 - t1:.3f}s")


# slider positions: the sorted union of metadata values from both datasets.
# At each position, show the dataset whose metadata[slider_key] is the largest
# value that is <= slider value (i.e. step-down / left-continuous lookup).
slider_positions = sorted(set(keys_left).union(keys_right))


def title_for(v, i_left, i_right):
    return (
        f"slider {slider_key} = {v:.3g}  |  "
        f"raw {slider_key}={keys_left[i_left]:.3g}  |  "
        f"normalized {slider_key}={keys_right[i_right]:.3g}"
    )


steps = []
for v in slider_positions:
    i_left = max(bisect.bisect_right(keys_left, float(v)) - 1, 0)
    i_right = max(bisect.bisect_right(keys_right, float(v)) - 1, 0)
    visibility = (
        [j == i_left for j in range(n_left)]
        + [j == i_right for j in range(n_right)]
    )
    steps.append(
        dict(
            method="update",
            args=[
                {"visible": visibility},
                {"title.text": title_for(float(v), i_left, i_right)},
            ],
            label=f"{float(v):.3g}",
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
        ticklen=0,
        # hide the per-tick step labels (still used for currentvalue display)
        font={"size": 1, "color": "rgba(0,0,0,0)"},
        steps=steps,
    )
]
t3 = time.perf_counter()
print(f"[timing] build {len(steps)} slider steps: {t3 - t2:.3f}s")


# build figure with two scenes side-by-side
fig = go.Figure(data=traces)
xy_range = [-40, 40]
common_axes = dict(
    xaxis_title="x",
    yaxis_title="y",
    zaxis_title="P contribution",
    xaxis=dict(range=xy_range),
    yaxis=dict(range=xy_range),
)
fig.update_layout(
    sliders=sliders,
    scene=dict(
        domain={"x": [0.0, 0.48], "y": [0.0, 1.0]},
        # zaxis=dict(range=[z_min_left, z_max_left]),
        **common_axes,
    ),
    scene2=dict(
        domain={"x": [0.52, 1.0], "y": [0.0, 1.0]},
        # zaxis=dict(range=[z_min_right, z_max_right]),
        **common_axes,
    ),
    title=title_for(float(slider_positions[0]), 0, 0),
    annotations=[
        dict(
            text="raw P",
            x=0.24, y=1.0, xref="paper", yref="paper",
            showarrow=False, font=dict(size=14),
        ),
        dict(
            text="normalized P",
            x=0.76, y=1.0, xref="paper", yref="paper",
            showarrow=False, font=dict(size=14),
        ),
    ],
)
t4 = time.perf_counter()
print(f"[timing] create figure + update_layout: {t4 - t3:.3f}s")


# Lock cameras between the two 3D scenes: when one rotates/zooms/pans, the
# other follows. Implemented by listening to plotly_relayout in the browser
# and mirroring the camera between scene <-> scene2.
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
