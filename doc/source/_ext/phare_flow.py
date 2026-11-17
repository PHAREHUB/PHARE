"""``.. flow::`` directive: a vertical chain of clickable boxes (plain HTML/CSS, no JavaScript, no image).

    .. flow::

       Python script | pyphare.pharein | /dev/python_to_cpp
       fills
       PHAREDict | src/initializer | /dev/python_to_cpp
       calls
       Simulator | src/simulator

Lines:   ``title | small text | target``  is a box (target optional: a document path such as ``/dev/x``, or ``https://...``)
         a line that does not contain ``|`` is the label of the arrow between two boxes (use ``-`` for an arrow without label)
         a blank line starts a new group; groups sit side by side on wide screens (one under the other on phones).
Option:  ``:horizontal:`` lay the boxes of a group out left to right (wraps on phones).
"""

from docutils import nodes
from docutils.parsers.rst import Directive, directives
from sphinx import addnodes


class Flow(Directive):
    has_content = True
    option_spec = {"horizontal": directives.flag}

    def _box(self, title, small, target):
        box = nodes.container(classes=["flow-box"] + (["flow-link"] if target else []))
        head = nodes.inline(text=title, classes=["flow-title"])
        if target:
            if target.startswith(("http://", "https://")):
                ref = nodes.reference("", "", refuri=target, internal=False)
                ref += head
            else:
                ref = addnodes.pending_xref("", head, refdomain="std", reftype="doc", reftarget=target,
                                            refexplicit=True, refwarn=True)
            box += nodes.paragraph("", "", ref)
        else:
            box += nodes.paragraph("", "", head)
        if small:
            box += nodes.paragraph(text=small, classes=["flow-small"])
        return box

    def run(self):
        groups, cur = [], []
        for raw in list(self.content) + [""]:
            line = raw.strip()
            if not line:
                if cur:
                    groups.append(cur)
                    cur = []
                continue
            cur.append(line)
        root = nodes.container(classes=["flow"] + (["flow-horizontal"] if "horizontal" in self.options else []))
        for g in groups:
            col = nodes.container(classes=["flow-group"])
            for line in g:
                if "|" in line:
                    parts = [p.strip() for p in line.split("|")]
                    parts += [""] * (3 - len(parts))
                    if len(col) and "flow-arrow" not in col[-1]["classes"]:
                        col += nodes.container(classes=["flow-arrow"])  # unlabeled arrow between consecutive boxes
                    col += self._box(*parts[:3])
                else:
                    label = "" if line == "-" else line
                    arrow = nodes.container(classes=["flow-arrow"])
                    arrow += nodes.paragraph(text=label)
                    col += arrow
            root += col
        return [root]


def setup(app):
    app.add_directive("flow", Flow)
    return {"parallel_read_safe": True, "version": "1"}
