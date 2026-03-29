"""
Privacy leakage audits through attack simulations.

Implements:
- MembershipInferenceAttack: Determine if a patient is in the database
- AttributeInferenceAttack: Infer hidden phenotypes from query results
- LeakageAuditReport: Comprehensive privacy evaluation
"""

import numpy as np
import pandas as pd
from typing import Dict, List, Set, Tuple, Optional, Any
from collections import defaultdict
import logging
from dataclasses import dataclass

try:
    from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
    from sklearn.model_selection import train_test_split
    from sklearn.metrics import (
        accuracy_score, precision_score, recall_score,
        f1_score, roc_auc_score, confusion_matrix
    )
    SKLEARN_AVAILABLE = True
except ImportError:
    SKLEARN_AVAILABLE = False

logger = logging.getLogger(__name__)


@dataclass
class AttackResults:
    """Container for attack evaluation results."""
    accuracy: float
    precision: float
    recall: float
    f1: float
    auc: float
    advantage: float  # Improvement over random guess
    confusion_matrix: Optional[np.ndarray] = None


class MembershipInferenceAttack:
    """
    Simulate membership inference attack.

    Goal: Given a target patient's phenotypes and access to the matching
    system, determine if the patient is in the database.

    Attack strategy:
    1. Train shadow models on known in/out data
    2. Extract features from query responses
    3. Train attack classifier on shadow data
    4. Evaluate on target database
    """

    def __init__(
        self,
        calculator,
        phenopackets: List[Dict],
        n_shadow_models: int = 5,
        shadow_fraction: float = 0.5
    ):
        """
        Initialize membership inference attack.

        Args:
            calculator: Similarity calculator to attack
            phenopackets: Full database of phenopackets
            n_shadow_models: Number of shadow models to train
            shadow_fraction: Fraction of data for each shadow model
        """
        if not SKLEARN_AVAILABLE:
            raise ImportError("scikit-learn required for attack simulations")

        self.calculator = calculator
        self.phenopackets = phenopackets
        self.n_shadow_models = n_shadow_models
        self.shadow_fraction = shadow_fraction

        self.attack_model = None
        self._is_trained = False

    def _extract_features(
        self,
        target: Dict,
        database: List[Dict],
        top_k: int = 10
    ) -> np.ndarray:
        """
        Extract features from query response for attack model.

        Features capture the distribution of similarity scores,
        which may differ for members vs non-members.

        Args:
            target: Target phenopacket to query
            database: Database to query against
            top_k: Number of top results to consider

        Returns:
            Feature vector
        """
        # Query the database
        matches = self.calculator.find_most_similar(target, database, top_k=top_k)

        # Handle suppressed results
        if matches is None:
            return np.zeros(10)  # Zero features if suppressed

        similarities = [sim for _, sim in matches]

        if not similarities:
            return np.zeros(10)

        # Extract statistical features
        features = [
            np.max(similarities),           # Max similarity
            np.mean(similarities),          # Mean similarity
            np.min(similarities),           # Min similarity
            np.std(similarities),           # Std similarity
            np.median(similarities),        # Median similarity
            similarities[0] - similarities[-1] if len(similarities) > 1 else 0,  # Gap
            len([s for s in similarities if s > 0.5]),  # High similarity count
            len([s for s in similarities if s > 0.3]),  # Medium similarity count
            np.percentile(similarities, 75) if len(similarities) >= 4 else np.max(similarities),  # 75th percentile
            np.percentile(similarities, 25) if len(similarities) >= 4 else np.min(similarities),  # 25th percentile
        ]

        return np.array(features)

    def train_attack_model(
        self,
        n_samples_per_class: int = 200,
        verbose: bool = False
    ):
        """
        Train attack model using shadow training.

        1. Create shadow databases (random subsets)
        2. Query with in/out samples
        3. Train classifier on (features, in_or_out)

        Args:
            n_samples_per_class: Samples per in/out class per shadow model
            verbose: Print progress
        """
        X_train = []
        y_train = []

        n_total = len(self.phenopackets)
        samples_per_shadow = n_samples_per_class // self.n_shadow_models

        for shadow_idx in range(self.n_shadow_models):
            if verbose:
                logger.info(f"Training shadow model {shadow_idx + 1}/{self.n_shadow_models}")

            # Create shadow database (random subset)
            shadow_size = int(n_total * self.shadow_fraction)
            shadow_indices = set(np.random.choice(n_total, shadow_size, replace=False))
            shadow_db = [self.phenopackets[i] for i in shadow_indices]

            # Sample "in" members
            in_indices = list(shadow_indices)[:samples_per_shadow]
            for idx in in_indices:
                target = self.phenopackets[idx]
                features = self._extract_features(target, shadow_db)
                X_train.append(features)
                y_train.append(1)  # In the database

            # Sample "out" non-members
            out_indices = [i for i in range(n_total) if i not in shadow_indices][:samples_per_shadow]
            for idx in out_indices:
                target = self.phenopackets[idx]
                features = self._extract_features(target, shadow_db)
                X_train.append(features)
                y_train.append(0)  # Not in the database

        X_train = np.array(X_train)
        y_train = np.array(y_train)

        # Train attack classifier
        self.attack_model = RandomForestClassifier(
            n_estimators=100,
            max_depth=10,
            random_state=42
        )
        self.attack_model.fit(X_train, y_train)
        self._is_trained = True

        if verbose:
            # Cross-validation estimate
            train_acc = self.attack_model.score(X_train, y_train)
            logger.info(f"Attack model training accuracy: {train_acc:.4f}")

    def evaluate_attack(
        self,
        n_test_samples: int = 200,
        test_fraction: float = 0.25
    ) -> AttackResults:
        """
        Evaluate attack success rate.

        Creates a test database and evaluates attack model's
        ability to distinguish members from non-members.

        Args:
            n_test_samples: Number of test samples per class
            test_fraction: Fraction of data for test database

        Returns:
            AttackResults with evaluation metrics
        """
        if not self._is_trained:
            raise ValueError("Attack model not trained. Call train_attack_model first.")

        n_total = len(self.phenopackets)

        # Create test database
        test_size = int(n_total * test_fraction)
        test_indices = set(np.random.choice(n_total, test_size, replace=False))
        test_db = [self.phenopackets[i] for i in test_indices]

        X_test = []
        y_test = []

        # Test on "in" samples
        in_indices = list(test_indices)[:n_test_samples]
        for idx in in_indices:
            target = self.phenopackets[idx]
            features = self._extract_features(target, test_db)
            X_test.append(features)
            y_test.append(1)

        # Test on "out" samples
        out_indices = [i for i in range(n_total) if i not in test_indices][:n_test_samples]
        for idx in out_indices:
            target = self.phenopackets[idx]
            features = self._extract_features(target, test_db)
            X_test.append(features)
            y_test.append(0)

        X_test = np.array(X_test)
        y_test = np.array(y_test)

        # Predict
        y_pred = self.attack_model.predict(X_test)
        y_prob = self.attack_model.predict_proba(X_test)[:, 1]

        # Compute metrics
        accuracy = accuracy_score(y_test, y_pred)
        precision = precision_score(y_test, y_pred, zero_division=0)
        recall = recall_score(y_test, y_pred, zero_division=0)
        f1 = f1_score(y_test, y_pred, zero_division=0)
        auc = roc_auc_score(y_test, y_prob)
        cm = confusion_matrix(y_test, y_pred)

        return AttackResults(
            accuracy=accuracy,
            precision=precision,
            recall=recall,
            f1=f1,
            auc=auc,
            advantage=accuracy - 0.5,  # Improvement over random guess
            confusion_matrix=cm
        )


class AttributeInferenceAttack:
    """
    Simulate attribute inference attack.

    Goal: Given partial phenotypes and query results, infer
    hidden/rare phenotypes that weren't directly queried.

    Attack strategy:
    1. Create censored queries (remove target term)
    2. Query database and get similar patients
    3. Check if similar patients have the target term
    4. Infer presence based on correlation
    """

    def __init__(
        self,
        calculator,
        phenopackets: List[Dict]
    ):
        """
        Initialize attribute inference attack.

        Args:
            calculator: Similarity calculator to attack
            phenopackets: Database of phenopackets
        """
        self.calculator = calculator
        self.phenopackets = phenopackets

        # Build term index for efficiency
        self._term_index = self._build_term_index()

    def _build_term_index(self) -> Dict[str, Set[int]]:
        """Build index of which patients have which terms."""
        index = defaultdict(set)
        for i, pp in enumerate(self.phenopackets):
            for feature in pp.get("phenotypicFeatures", []):
                if not feature.get("excluded", False):
                    term_id = feature["type"]["id"]
                    index[term_id].add(i)
        return index

    def _get_patient_terms(self, phenopacket: Dict) -> Set[str]:
        """Get set of HPO terms from phenopacket."""
        terms = set()
        for feature in phenopacket.get("phenotypicFeatures", []):
            if not feature.get("excluded", False):
                terms.add(feature["type"]["id"])
        return terms

    def _censor_term(self, phenopacket: Dict, term_id: str) -> Dict:
        """Create copy of phenopacket with term removed."""
        import copy
        censored = copy.deepcopy(phenopacket)
        censored["phenotypicFeatures"] = [
            f for f in censored.get("phenotypicFeatures", [])
            if f["type"]["id"] != term_id
        ]
        return censored

    def attack_term(
        self,
        target_term: str,
        n_trials: int = 100,
        top_k: int = 5,
        inference_threshold: int = 3
    ) -> Dict[str, Any]:
        """
        Attempt to infer presence of a specific term.

        Strategy: If top-k similar patients mostly have the term,
        predict the target also has it.

        Args:
            target_term: HPO term ID to infer
            n_trials: Number of patients to test
            top_k: Number of similar patients to check
            inference_threshold: Number of matches needed to predict positive

        Returns:
            Attack evaluation metrics
        """
        # Get patients with this term
        patients_with_term = self._term_index.get(target_term, set())

        if len(patients_with_term) < 10:
            logger.warning(f"Insufficient patients with term {target_term}")
            return {"error": "insufficient_samples"}

        y_true = []
        y_pred = []

        # Test on random sample of patients
        test_indices = np.random.choice(
            len(self.phenopackets),
            min(n_trials, len(self.phenopackets)),
            replace=False
        )

        for idx in test_indices:
            pp = self.phenopackets[idx]
            patient_terms = self._get_patient_terms(pp)

            # Ground truth: does patient have the term?
            has_term = target_term in patient_terms
            y_true.append(int(has_term))

            # Create censored query (without target term)
            censored = self._censor_term(pp, target_term)

            # Query for similar patients
            matches = self.calculator.find_most_similar(
                censored, self.phenopackets, top_k=top_k
            )

            if matches is None:
                # Suppressed - predict negative
                y_pred.append(0)
                continue

            # Check how many similar patients have the term
            term_count = 0
            for match_idx, _ in matches:
                if match_idx == idx:
                    continue  # Skip self
                match_terms = self._get_patient_terms(self.phenopackets[match_idx])
                if target_term in match_terms:
                    term_count += 1

            # Inference: if many similar patients have term, predict positive
            predicted = int(term_count >= inference_threshold)
            y_pred.append(predicted)

        # Compute metrics
        y_true = np.array(y_true)
        y_pred = np.array(y_pred)

        if len(np.unique(y_true)) < 2:
            return {"error": "no_variance_in_labels"}

        return {
            "term": target_term,
            "n_patients_with_term": len(patients_with_term),
            "n_trials": len(y_true),
            "baseline_rate": y_true.mean(),
            "accuracy": accuracy_score(y_true, y_pred),
            "precision": precision_score(y_true, y_pred, zero_division=0),
            "recall": recall_score(y_true, y_pred, zero_division=0),
            "f1": f1_score(y_true, y_pred, zero_division=0),
            "advantage": accuracy_score(y_true, y_pred) - max(y_true.mean(), 1 - y_true.mean())
        }

    def attack_rare_terms(
        self,
        max_prevalence: float = 0.2,
        min_prevalence: float = 0.05,
        n_terms: int = 5,
        n_trials_per_term: int = 50
    ) -> pd.DataFrame:
        """
        Attack multiple rare terms and aggregate results.

        Args:
            max_prevalence: Maximum term prevalence to consider "rare"
            min_prevalence: Minimum prevalence (need enough samples)
            n_terms: Number of terms to attack
            n_trials_per_term: Trials per term

        Returns:
            DataFrame with attack results per term
        """
        n_patients = len(self.phenopackets)

        # Find rare terms within prevalence range
        rare_terms = []
        for term, indices in self._term_index.items():
            prevalence = len(indices) / n_patients
            if min_prevalence <= prevalence <= max_prevalence:
                rare_terms.append((term, prevalence))

        # Sort by prevalence and take subset
        rare_terms.sort(key=lambda x: x[1])
        selected_terms = rare_terms[:n_terms]

        results = []
        for term, prevalence in selected_terms:
            logger.info(f"Attacking term {term} (prevalence={prevalence:.3f})")
            result = self.attack_term(term, n_trials=n_trials_per_term)
            if "error" not in result:
                result["prevalence"] = prevalence
                results.append(result)

        return pd.DataFrame(results)


class LeakageAuditReport:
    """
    Generate comprehensive privacy leakage audit report.

    Compares attack success rates between:
    - Baseline (no privacy) system
    - Privacy-preserving system

    Quantifies privacy gain from each mechanism.
    """

    def __init__(
        self,
        base_calculator,
        private_calculator,
        phenopackets: List[Dict]
    ):
        """
        Initialize leakage audit.

        Args:
            base_calculator: Non-private calculator (baseline)
            private_calculator: Privacy-preserving calculator
            phenopackets: Database to audit
        """
        self.base_calculator = base_calculator
        self.private_calculator = private_calculator
        self.phenopackets = phenopackets

    def run_membership_inference_audit(
        self,
        n_shadow_models: int = 5,
        n_test_samples: int = 100,
        verbose: bool = True
    ) -> Dict[str, Any]:
        """
        Audit membership inference vulnerability.

        Args:
            n_shadow_models: Shadow models for training
            n_test_samples: Test samples per class
            verbose: Print progress

        Returns:
            Audit results comparing baseline and private
        """
        results = {}

        # Baseline attack
        if verbose:
            logger.info("Running membership inference on baseline system...")

        mia_base = MembershipInferenceAttack(
            self.base_calculator,
            self.phenopackets,
            n_shadow_models=n_shadow_models
        )
        mia_base.train_attack_model(verbose=verbose)
        base_results = mia_base.evaluate_attack(n_test_samples=n_test_samples)

        results["baseline"] = {
            "accuracy": base_results.accuracy,
            "precision": base_results.precision,
            "recall": base_results.recall,
            "f1": base_results.f1,
            "auc": base_results.auc,
            "advantage": base_results.advantage
        }

        # Private system attack
        if verbose:
            logger.info("Running membership inference on private system...")

        mia_private = MembershipInferenceAttack(
            self.private_calculator,
            self.phenopackets,
            n_shadow_models=n_shadow_models
        )
        mia_private.train_attack_model(verbose=verbose)
        private_results = mia_private.evaluate_attack(n_test_samples=n_test_samples)

        results["private"] = {
            "accuracy": private_results.accuracy,
            "precision": private_results.precision,
            "recall": private_results.recall,
            "f1": private_results.f1,
            "auc": private_results.auc,
            "advantage": private_results.advantage
        }

        # Privacy gain
        results["privacy_gain"] = {
            "accuracy_reduction": base_results.accuracy - private_results.accuracy,
            "advantage_reduction": base_results.advantage - private_results.advantage,
            "auc_reduction": base_results.auc - private_results.auc
        }

        return results

    def run_attribute_inference_audit(
        self,
        n_terms: int = 5,
        n_trials_per_term: int = 50,
        verbose: bool = True
    ) -> Dict[str, Any]:
        """
        Audit attribute inference vulnerability.

        Args:
            n_terms: Number of rare terms to attack
            n_trials_per_term: Trials per term
            verbose: Print progress

        Returns:
            Audit results comparing baseline and private
        """
        results = {}

        # Baseline attack
        if verbose:
            logger.info("Running attribute inference on baseline system...")

        aia_base = AttributeInferenceAttack(
            self.base_calculator,
            self.phenopackets
        )
        base_df = aia_base.attack_rare_terms(
            n_terms=n_terms,
            n_trials_per_term=n_trials_per_term
        )

        if len(base_df) > 0:
            results["baseline"] = {
                "mean_accuracy": base_df["accuracy"].mean(),
                "mean_advantage": base_df["advantage"].mean(),
                "max_advantage": base_df["advantage"].max(),
                "per_term": base_df.to_dict("records")
            }
        else:
            results["baseline"] = {"error": "no_valid_terms"}

        # Private system attack
        if verbose:
            logger.info("Running attribute inference on private system...")

        aia_private = AttributeInferenceAttack(
            self.private_calculator,
            self.phenopackets
        )
        private_df = aia_private.attack_rare_terms(
            n_terms=n_terms,
            n_trials_per_term=n_trials_per_term
        )

        if len(private_df) > 0:
            results["private"] = {
                "mean_accuracy": private_df["accuracy"].mean(),
                "mean_advantage": private_df["advantage"].mean(),
                "max_advantage": private_df["advantage"].max(),
                "per_term": private_df.to_dict("records")
            }
        else:
            results["private"] = {"error": "no_valid_terms"}

        # Privacy gain
        if "error" not in results["baseline"] and "error" not in results["private"]:
            results["privacy_gain"] = {
                "mean_advantage_reduction": (
                    results["baseline"]["mean_advantage"] -
                    results["private"]["mean_advantage"]
                ),
                "max_advantage_reduction": (
                    results["baseline"]["max_advantage"] -
                    results["private"]["max_advantage"]
                )
            }

        return results

    def run_full_audit(
        self,
        verbose: bool = True
    ) -> Dict[str, Any]:
        """
        Run comprehensive privacy audit.

        Args:
            verbose: Print progress

        Returns:
            Complete audit results
        """
        if verbose:
            logger.info("Starting comprehensive privacy audit...")

        results = {
            "n_phenopackets": len(self.phenopackets)
        }

        # Membership inference
        results["membership_inference"] = self.run_membership_inference_audit(
            verbose=verbose
        )

        # Attribute inference
        results["attribute_inference"] = self.run_attribute_inference_audit(
            verbose=verbose
        )

        # Summary
        mia_gain = results["membership_inference"].get("privacy_gain", {})
        aia_gain = results["attribute_inference"].get("privacy_gain", {})

        results["summary"] = {
            "mia_advantage_reduction": mia_gain.get("advantage_reduction", 0),
            "mia_auc_reduction": mia_gain.get("auc_reduction", 0),
            "aia_advantage_reduction": aia_gain.get("mean_advantage_reduction", 0),
            "privacy_improved": (
                mia_gain.get("advantage_reduction", 0) > 0 or
                aia_gain.get("mean_advantage_reduction", 0) > 0
            )
        }

        if verbose:
            logger.info("Privacy audit complete.")
            logger.info(f"MIA advantage reduction: {results['summary']['mia_advantage_reduction']:.4f}")
            logger.info(f"AIA advantage reduction: {results['summary']['aia_advantage_reduction']:.4f}")

        return results


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)

    if SKLEARN_AVAILABLE:
        print("Leakage audit module loaded successfully")
        print("scikit-learn available for attack simulations")
    else:
        print("WARNING: scikit-learn not available")
        print("Install with: pip install scikit-learn")
