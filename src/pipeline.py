"""Publication Bias Forensics Pipeline.

Runs all 8 bias methods on 403 Cochrane reviews from the Fragility Atlas dataset.
Usage: python -m src.pipeline [--pairwise-dir DIR] [--output-dir DIR]
"""

import sys
import json
import time
import csv
import argparse
import numpy as np
from pathlib import Path

from src.methods import run_all_methods
from src.loader import load_all_reviews


DEFAULT_PAIRWISE_DIR = r'C:\Models\Pairwise70\data'
DEFAULT_OUTPUT_DIR = r'C:\BiasForensics\data\output'


def classify_bias(fingerprint):
    """Classify a review's bias profile as Clean/Suspected/Confirmed/Discordant."""
    # Count detection tests that are significant
    n_detect = 0
    if fingerprint['egger']['significant']:
        n_detect += 1
    if fingerprint['begg']['significant']:
        n_detect += 1
    if fingerprint['pcurve']['inadequate']:
        n_detect += 1
    if fingerprint['sel3psm']['significant']:
        n_detect += 1

    # Check if correction methods substantially change the effect
    theta_unadj = fingerprint['trimfill']['theta_unadj']
    corrections = [
        fingerprint['trimfill']['theta_adj'],
        fingerprint['petpeese']['theta_adj'],
        fingerprint['sel3psm']['theta_adj'],
        fingerprint['limit']['theta_limit'],
    ]

    # Relative shift: how much do corrections change the effect?
    shifts = []
    for c in corrections:
        if abs(theta_unadj) > 0.01:
            shift = abs(c - theta_unadj) / abs(theta_unadj)
        else:
            shift = abs(c - theta_unadj)
        shifts.append(shift)

    max_shift = max(shifts) if shifts else 0
    mean_shift = np.mean(shifts) if shifts else 0

    # Do corrections agree on direction?
    directions = [1 if c >= 0 else -1 for c in corrections]
    unadj_dir = 1 if theta_unadj >= 0 else -1
    same_dir = all(d == unadj_dir for d in directions)

    # Do corrections agree on significance?
    # (simplified: just check if they're all on the same side of zero)
    n_agree_dir = sum(1 for d in directions if d == unadj_dir)

    # Classification
    if n_detect <= 1 and mean_shift < 0.2:
        bias_class = 'Clean'
    elif n_detect >= 3 and mean_shift >= 0.2:
        bias_class = 'Confirmed'
    elif n_agree_dir <= 2:  # corrections point in different directions
        bias_class = 'Discordant'
    else:
        bias_class = 'Suspected'

    return {
        'n_detect': n_detect,
        'max_shift': round(max_shift, 4),
        'mean_shift': round(mean_shift, 4),
        'concordance': round(n_agree_dir / 4 * 100, 1),
        'bias_class': bias_class,
    }


def run_pipeline(pairwise_dir, output_dir, max_reviews=0):
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    print("Publication Bias Forensics Pipeline")
    print("=" * 40)
    print(f"Data: {pairwise_dir}")
    print()

    # Load reviews
    print("Loading reviews...")
    reviews = list(load_all_reviews(pairwise_dir, min_k=3))
    # Filter to k >= 5 for meaningful bias tests
    reviews = [r for r in reviews if r.k >= 5]
    if max_reviews > 0:
        reviews = reviews[:max_reviews]
    print(f"  {len(reviews)} reviews with k >= 5")
    print()

    # Run forensics
    print("Running 8 bias methods on each review...")
    results = []
    t0 = time.time()

    for i, review in enumerate(reviews):
        try:
            fingerprint = run_all_methods(review.yi, review.sei)
            classification = classify_bias(fingerprint)

            row = {
                'review_id': review.review_id,
                'k': review.k,
                'scale': review.scale,
                'analysis_name': review.analysis_name,
                # Egger
                'egger_p': round(fingerprint['egger']['p'], 4),
                'egger_sig': int(fingerprint['egger']['significant']),
                # Begg
                'begg_p': round(fingerprint['begg']['p'], 4),
                'begg_tau': round(fingerprint['begg']['tau'], 4),
                'begg_sig': int(fingerprint['begg']['significant']),
                # P-curve
                'pcurve_n_sig': fingerprint['pcurve']['n_sig'],
                'pcurve_evidential': int(fingerprint['pcurve']['evidential']),
                'pcurve_inadequate': int(fingerprint['pcurve']['inadequate']),
                # Z-curve
                'zcurve_oir': round(fingerprint['zcurve']['oir'], 4),
                'zcurve_edr': round(fingerprint['zcurve']['edr'], 4),
                # Trim-and-fill
                'tf_theta_unadj': round(fingerprint['trimfill']['theta_unadj'], 4),
                'tf_theta_adj': round(fingerprint['trimfill']['theta_adj'], 4),
                'tf_k0': fingerprint['trimfill']['k0'],
                # PET-PEESE
                'petpeese_theta': round(fingerprint['petpeese']['theta_adj'], 4),
                'petpeese_method': fingerprint['petpeese']['method_used'],
                # Selection model
                'sel3psm_theta': round(fingerprint['sel3psm']['theta_adj'], 4),
                'sel3psm_eta': round(fingerprint['sel3psm']['eta'], 4),
                'sel3psm_lr_p': round(fingerprint['sel3psm']['lr_p'], 4),
                'sel3psm_sig': int(fingerprint['sel3psm']['significant']),
                # Limit meta
                'limit_theta': round(fingerprint['limit']['theta_limit'], 4),
                # Classification
                'n_detect': classification['n_detect'],
                'max_shift': classification['max_shift'],
                'mean_shift': classification['mean_shift'],
                'concordance': classification['concordance'],
                'bias_class': classification['bias_class'],
            }
            results.append(row)
        except Exception as e:
            print(f"  Warning: {review.review_id} failed: {e}", file=sys.stderr)
            continue

        if (i + 1) % 20 == 0 or (i + 1) == len(reviews):
            elapsed = time.time() - t0
            rate = (i + 1) / elapsed if elapsed > 0 else 0
            eta = (len(reviews) - i - 1) / rate if rate > 0 else 0
            print(f"  [{i+1}/{len(reviews)}] {review.review_id} ({review.k} studies) "
                  f"-> {row['bias_class']} [{elapsed:.0f}s, ETA {eta:.0f}s]")
            sys.stdout.flush()

    total_time = time.time() - t0
    print(f"\n  Done: {len(results)} reviews in {total_time:.1f}s")
    print()

    # Export
    results_path = output_path / 'bias_forensics_results.csv'
    if results:
        fields = list(results[0].keys())
        with open(results_path, 'w', newline='', encoding='utf-8') as f:
            writer = csv.DictWriter(f, fieldnames=fields)
            writer.writeheader()
            for r in results:
                writer.writerow(r)
    print(f"  Results: {results_path}")

    # Summary
    n = len(results)
    if n > 0:
        counts = {}
        for cat in ['Clean', 'Suspected', 'Confirmed', 'Discordant']:
            counts[cat] = sum(1 for r in results if r['bias_class'] == cat)

        summary = {
            'n_reviews': n,
            'classification_counts': counts,
            'method_detection_rates': {
                'egger': round(sum(r['egger_sig'] for r in results) / n * 100, 1),
                'begg': round(sum(r['begg_sig'] for r in results) / n * 100, 1),
                'pcurve_inadequate': round(sum(r['pcurve_inadequate'] for r in results) / n * 100, 1),
                'selection_model': round(sum(r['sel3psm_sig'] for r in results) / n * 100, 1),
            },
            'mean_tf_k0': round(np.mean([r['tf_k0'] for r in results]), 2),
            'mean_shift': round(np.mean([r['mean_shift'] for r in results]), 4),
            'elapsed_seconds': round(total_time, 1),
        }

        summary_path = output_path / 'bias_forensics_summary.json'
        with open(summary_path, 'w', encoding='utf-8') as f:
            json.dump(summary, f, indent=2)
        print(f"  Summary: {summary_path}")

        # Headline
        print()
        print("=" * 50)
        print("HEADLINE RESULTS")
        print("=" * 50)
        for cat in ['Clean', 'Suspected', 'Confirmed', 'Discordant']:
            pct = counts[cat] / n * 100
            print(f"  {cat:12s}: {counts[cat]:4d} ({pct:5.1f}%)")
        print()
        print(f"  Detection rates:")
        for method, rate in summary['method_detection_rates'].items():
            print(f"    {method:20s}: {rate:5.1f}%")
        print(f"  Mean trim-and-fill k0: {summary['mean_tf_k0']:.1f}")
        print(f"  Mean effect shift: {summary['mean_shift']:.3f}")

    return results


def main():
    parser = argparse.ArgumentParser(description='Publication Bias Forensics')
    parser.add_argument('--pairwise-dir', default=DEFAULT_PAIRWISE_DIR)
    parser.add_argument('--output-dir', default=DEFAULT_OUTPUT_DIR)
    parser.add_argument('--max-reviews', type=int, default=0)
    args = parser.parse_args()
    run_pipeline(args.pairwise_dir, args.output_dir, args.max_reviews)


if __name__ == '__main__':
    main()
