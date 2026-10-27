"""Check plotted coordinates, interval geometry, quartiles, and page layout."""
import csv
import os
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parent.parent
os.environ.setdefault('MPLCONFIGDIR', str(ROOT / '.mplconfig'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.dates as mdates
import numpy as np
import pandas as pd
import config
from sarflow.temporal import validate
from sarflow.publication import plot_three_pair_coherence


class PublicationTests(unittest.TestCase):
    def test_actual_statistics_geometry_and_layout(self):
        with (Path(config.Output_Folder)/'timeseries_summary.csv').open(newline='') as stream:
            records = list(csv.DictReader(stream))
        table = pd.DataFrame([{key: value if key == 'observation' else float(value)
                               for key, value in row.items()} for row in records])
        obs = validate(config)
        dates = mdates.date2num([o.end for o in obs])

        def inspect(fig, out, stem):
            ax = fig.axes[0]
            self.assertEqual(ax.get_title(), '')
            self.assertEqual(ax.get_xlabel(), 'Pair end date (UTC)')
            np.testing.assert_array_equal(ax.lines[0].get_xdata(), dates)
            np.testing.assert_array_equal(ax.lines[0].get_ydata(), table.mean_exact)
            np.testing.assert_array_equal(ax.lines[1].get_xdata(), dates)
            np.testing.assert_array_equal(ax.lines[1].get_ydata(), table.median_sample)
            self.assertEqual([ax.lines[i].get_marker() for i in (2, 6, 10)], ['o', 's', '^'])
            band = ax.collections[0].get_paths()[0].vertices
            for i, date in enumerate(dates):
                vals = band[np.isclose(band[:, 0], date), 1]
                self.assertTrue(np.any(vals == table.q25_sample.iloc[i]))
                self.assertTrue(np.any(vals == table.q75_sample.iloc[i]))
            for container, o in zip(ax.containers, obs):
                segment = container.lines[2][0].get_segments()[0]
                np.testing.assert_array_equal(segment[:, 0], mdates.date2num([o.start, o.end]))
            np.testing.assert_array_equal(ax.lines[-1].get_xdata(), [mdates.date2num(obs[0].end)+6]*2)
            tick_days = np.diff(ax.get_xticks())
            self.assertTrue((tick_days == tick_days[0]).all())
            for label in ax.get_xticklabels():
                self.assertRegex(label.get_text(), r'^\d{4}-\d{2}-\d{2}$')
            self.assertFalse(ax.spines['top'].get_visible())
            self.assertFalse(ax.spines['right'].get_visible())
            fig.canvas.draw()
            renderer = fig.canvas.get_renderer()
            boxes = [label.get_window_extent(renderer) for label in ax.get_xticklabels()]
            for left, right in zip(boxes, boxes[1:]):
                self.assertLess(left.x1+2, right.x0)
            artists = [ax.xaxis.label, ax.yaxis.label, ax.get_legend(), *ax.texts, *fig.texts,
                       *ax.get_xticklabels(), *ax.get_yticklabels()]
            page = fig.bbox
            for artist in artists:
                box = artist.get_window_extent(renderer)
                self.assertGreaterEqual(box.x0, page.x0-.5)
                self.assertGreaterEqual(box.y0, page.y0-.5)
                self.assertLessEqual(box.x1, page.x1+.5)
                self.assertLessEqual(box.y1, page.y1+.5)

        with TemporaryDirectory(dir=ROOT) as folder:
            with patch('sarflow.publication.export_publication_figure', side_effect=inspect):
                plot_three_pair_coherence(table, obs, Path(folder), config)
                plot_three_pair_coherence(table, obs, Path(folder), config, single_column=True)


if __name__ == '__main__':
    unittest.main(verbosity=2)
