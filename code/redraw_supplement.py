"""Export the five checked SI figures without a manuscript dependency."""
import hashlib
import json
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.text import Text
from PIL import Image
import plot_supplement_functions as si


def main():
    si.setup()
    source = si.TEX.read_text(encoding="utf-8")
    spec = [
        ("supplement_1_lifecycle", *si.si1(json.loads(si.STRUCTURAL.read_text()))),
        ("supplement_2_ownership", *si.si2(si.read_csv(si.OWNERSHIP))),
        ("supplement_3_collection", *si.si3(si.read_csv(si.COLLECTION))),
        ("supplement_4_events", *si.si4(source)),
        ("supplement_5_attainment", *si.si5(source)),
    ]
    manifest = dict(native_raster_dpi=1200, vectors_have_no_dpi=True, figures=[])
    si.OUT.mkdir(exist_ok=True, parents=True)
    for name, fig, data in spec:
        audit = si.layout_audit(fig)
        assert audit["zero_detected_overlap_or_clipping"], (name, audit)
        points = [x.get_fontsize() for x in fig.findobj(match=Text) if x.get_visible() and x.get_text().strip()]
        assert min(points) * min(1, 164 / (fig.get_size_inches()[0] * 25.4)) >= 5
        outputs = []
        for ext in ("pdf", "svg", "png", "tiff"):
            target = si.OUT / f"{name}.{ext}"
            if target.exists():
                raise FileExistsError(target)
            options = dict(dpi=1200) if ext in ("png", "tiff") else {}
            if ext == "tiff":
                options["pil_kwargs"] = dict(compression="tiff_lzw")
            fig.savefig(target, facecolor="white", **options)
            if ext in ("png", "tiff"):
                with Image.open(target) as image:
                    assert all(abs(v - 1200) < .02 for v in image.info["dpi"])
            outputs.append(dict(file=target.name, sha256=hashlib.sha256(target.read_bytes()).hexdigest()))
        manifest["figures"].append(dict(name=name, data=data, layout_audit=audit, outputs=outputs))
        plt.close(fig)
        print(name + ": 4 exports; layout audit passed", flush=True)
    target = si.ROOT / "results/epse_v39/public_supplement_figure_manifest.json"
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
