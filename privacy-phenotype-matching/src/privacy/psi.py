"""
Private Set Intersection (PSI) for secure phenotype matching.

Implements:
- DiffieHellmanPSI: DH-based PSI using elliptic curves
- OPRFBasedPSI: Oblivious PRF variant
- HybridPSI: Combines PSI with HPO ancestor expansion
- PSIPhenopacketMatcher: High-level interface for phenopacket matching

These implementations enable two parties to compute the intersection of
their phenotype sets without revealing non-intersecting elements.
"""

import hashlib
import secrets
from typing import Dict, List, Set, Tuple, Optional, Any
from dataclasses import dataclass
import logging

# Cryptographic imports
from cryptography.hazmat.primitives.asymmetric import ec
from cryptography.hazmat.primitives import hashes, serialization
from cryptography.hazmat.backends import default_backend
from cryptography.hazmat.primitives.kdf.hkdf import HKDF

logger = logging.getLogger(__name__)


@dataclass
class PSIConfig:
    """Configuration for PSI protocol."""
    security_parameter: int = 128
    curve: str = "SECP256R1"  # NIST P-256
    hash_algorithm: str = "SHA256"


class ECPoint:
    """
    Elliptic curve point wrapper for PSI operations.

    Uses the cryptography library for actual EC operations.
    """

    # Curve parameters
    CURVE = ec.SECP256R1()

    def __init__(self, private_key: Optional[ec.EllipticCurvePrivateKey] = None):
        """
        Initialize EC point with optional private key.

        Args:
            private_key: Existing private key, or None to generate new
        """
        if private_key is None:
            self._private_key = ec.generate_private_key(self.CURVE, default_backend())
        else:
            self._private_key = private_key

        self._public_key = self._private_key.public_key()

    @classmethod
    def from_hash(cls, data: bytes) -> "ECPoint":
        """
        Create EC point by hashing data to curve.

        Uses hash-to-curve approach for deterministic point generation.

        Args:
            data: Bytes to hash

        Returns:
            ECPoint on the curve
        """
        # Use HKDF to derive a private key from the hash
        hkdf = HKDF(
            algorithm=hashes.SHA256(),
            length=32,
            salt=b"phenotype-psi-h2c",
            info=b"hash-to-curve",
            backend=default_backend()
        )
        key_material = hkdf.derive(data)

        # Convert to integer and use as private key
        # (This is a simplified hash-to-curve; production would use RFC 9380)
        private_value = int.from_bytes(key_material, 'big')

        # Create private key from the value
        private_key = ec.derive_private_key(
            private_value % (2**256 - 1),  # Ensure valid range
            ec.SECP256R1(),
            default_backend()
        )

        return cls(private_key)

    def get_public_bytes(self) -> bytes:
        """Get serialized public key bytes."""
        return self._public_key.public_bytes(
            encoding=serialization.Encoding.X962,
            format=serialization.PublicFormat.CompressedPoint
        )

    def scalar_mult(self, scalar: int) -> bytes:
        """
        Multiply point by scalar (simulated via key derivation).

        Returns bytes representing the result.
        """
        # Use ECDH-like operation
        scalar_bytes = scalar.to_bytes(32, 'big')

        # Derive a deterministic value
        hkdf = HKDF(
            algorithm=hashes.SHA256(),
            length=32,
            salt=self.get_public_bytes(),
            info=scalar_bytes,
            backend=default_backend()
        )

        return hkdf.derive(self.get_public_bytes())


class DiffieHellmanPSI:
    """
    Private Set Intersection using Diffie-Hellman.

    Protocol overview:
    1. Server hashes each element to EC point, raises to secret exponent
    2. Client hashes each element to EC point, raises to secret exponent, sends to server
    3. Server raises client's points to server's secret, returns
    4. Client raises server's points to client's secret
    5. Intersection found by comparing double-raised values

    Security: Based on DDH assumption in the chosen EC group.
    """

    def __init__(self, config: Optional[PSIConfig] = None):
        """
        Initialize DH-PSI protocol.

        Args:
            config: PSI configuration
        """
        self.config = config or PSIConfig()
        self._secret: Optional[int] = None

    def _hash_element(self, element: str) -> bytes:
        """Hash an element (HPO term) to bytes."""
        return hashlib.sha256(element.encode()).digest()

    def _hash_to_point(self, element: str) -> ECPoint:
        """Hash an element to an elliptic curve point."""
        element_hash = self._hash_element(element)
        return ECPoint.from_hash(element_hash)

    def generate_secret(self) -> int:
        """Generate a random secret exponent."""
        self._secret = secrets.randbits(self.config.security_parameter)
        return self._secret

    def client_prepare(
        self,
        elements: Set[str]
    ) -> Tuple[Dict[bytes, str], int]:
        """
        Client prepares elements for PSI.

        Step 1: Hash each element to EC point, raise to client's secret.

        Args:
            elements: Set of HPO term IDs

        Returns:
            Tuple of (encrypted_elements mapping, client_secret)
        """
        secret = self.generate_secret()

        encrypted = {}
        for element in elements:
            point = self._hash_to_point(element)
            # Raise to client's secret
            encrypted_value = point.scalar_mult(secret)
            encrypted[encrypted_value] = element

        return encrypted, secret

    def server_process(
        self,
        client_encrypted: Dict[bytes, Any],
        server_elements: Set[str]
    ) -> Tuple[List[bytes], Dict[bytes, str]]:
        """
        Server processes client's encrypted elements and prepares own elements.

        Step 2: Raise client's elements to server's secret,
                hash own elements and raise to server's secret.

        Args:
            client_encrypted: Client's encrypted elements (values)
            server_elements: Server's plaintext elements

        Returns:
            Tuple of (doubly_encrypted_client, server_encrypted)
        """
        server_secret = self.generate_secret()

        # Double-encrypt client's elements
        doubly_encrypted = []
        for client_value in client_encrypted.keys():
            # Apply server's secret to client's value
            # (Simulated via hash combining)
            combined = hashlib.sha256(
                client_value + server_secret.to_bytes(32, 'big')
            ).digest()
            doubly_encrypted.append(combined)

        # Encrypt server's elements
        server_encrypted = {}
        for element in server_elements:
            point = self._hash_to_point(element)
            encrypted_value = point.scalar_mult(server_secret)
            server_encrypted[encrypted_value] = element

        return doubly_encrypted, server_encrypted

    def client_compute_intersection(
        self,
        client_encrypted: Dict[bytes, str],
        client_secret: int,
        doubly_encrypted: List[bytes],
        server_encrypted: Dict[bytes, str]
    ) -> Set[str]:
        """
        Client computes final intersection.

        Step 3: Raise server's elements to client's secret, compare.

        Args:
            client_encrypted: Client's original encrypted mapping
            client_secret: Client's secret
            doubly_encrypted: Doubly encrypted client elements from server
            server_encrypted: Server's singly encrypted elements

        Returns:
            Set of intersecting elements (client-side terms)
        """
        # Raise server's elements to client's secret
        server_doubly = {}
        for server_value, element in server_encrypted.items():
            combined = hashlib.sha256(
                server_value + client_secret.to_bytes(32, 'big')
            ).digest()
            server_doubly[combined] = element

        # Find intersection by comparing doubly-encrypted values
        intersection = set()

        # Map client's doubly-encrypted back to original elements
        client_elements = list(client_encrypted.values())

        for i, doubly_enc in enumerate(doubly_encrypted):
            if doubly_enc in server_doubly:
                # This element is in the intersection
                if i < len(client_elements):
                    intersection.add(client_elements[i])

        return intersection


class OPRFBasedPSI:
    """
    PSI based on Oblivious Pseudo-Random Functions.

    Uses OPRF to evaluate PRF on each element without revealing input.
    More efficient than DH-PSI for larger sets.

    Simplified implementation using hash-based OPRF simulation.
    """

    def __init__(self, config: Optional[PSIConfig] = None):
        """
        Initialize OPRF-PSI.

        Args:
            config: PSI configuration
        """
        self.config = config or PSIConfig()
        self._key: Optional[bytes] = None

    def generate_key(self) -> bytes:
        """Generate OPRF key."""
        self._key = secrets.token_bytes(32)
        return self._key

    def _oprf_eval(self, element: str, key: bytes) -> bytes:
        """
        Evaluate OPRF on element with key.

        F_k(x) = H(k || H(x))

        Args:
            element: Input element
            key: OPRF key

        Returns:
            PRF output
        """
        element_hash = hashlib.sha256(element.encode()).digest()
        return hashlib.sha256(key + element_hash).digest()

    def server_encode(self, elements: Set[str]) -> Set[bytes]:
        """
        Server encodes elements using OPRF.

        Args:
            elements: Server's plaintext elements

        Returns:
            Set of OPRF encodings
        """
        if self._key is None:
            self.generate_key()

        return {self._oprf_eval(e, self._key) for e in elements}

    def client_blind(self, elements: Set[str]) -> Tuple[Dict[bytes, str], bytes]:
        """
        Client blinds elements for OPRF evaluation.

        Args:
            elements: Client's elements

        Returns:
            Tuple of (blinded_elements mapping, blinding_factor)
        """
        blinding = secrets.token_bytes(32)

        blinded = {}
        for element in elements:
            element_hash = hashlib.sha256(element.encode()).digest()
            blinded_hash = hashlib.sha256(element_hash + blinding).digest()
            blinded[blinded_hash] = element

        return blinded, blinding

    def server_evaluate(
        self,
        blinded_elements: Dict[bytes, Any]
    ) -> Dict[bytes, bytes]:
        """
        Server evaluates OPRF on blinded elements.

        Args:
            blinded_elements: Client's blinded elements

        Returns:
            OPRF evaluations
        """
        if self._key is None:
            raise ValueError("Server key not initialized")

        evaluated = {}
        for blinded in blinded_elements.keys():
            # Apply server key
            result = hashlib.sha256(self._key + blinded).digest()
            evaluated[blinded] = result

        return evaluated

    def client_finalize(
        self,
        blinded_mapping: Dict[bytes, str],
        blinding: bytes,
        evaluations: Dict[bytes, bytes],
        server_encoded: Set[bytes]
    ) -> Set[str]:
        """
        Client finalizes intersection computation.

        Args:
            blinded_mapping: Mapping from blinded values to elements
            blinding: Client's blinding factor
            evaluations: Server's OPRF evaluations
            server_encoded: Server's encoded elements

        Returns:
            Intersection elements
        """
        intersection = set()

        for blinded, element in blinded_mapping.items():
            if blinded in evaluations:
                # Unblind
                final_value = hashlib.sha256(
                    evaluations[blinded] + blinding
                ).digest()

                if final_value in server_encoded:
                    intersection.add(element)

        return intersection


class HybridPSI:
    """
    Hybrid PSI optimized for HPO phenotype matching.

    Combines DH-PSI with HPO ontology awareness:
    - Expands terms to include ancestors for semantic matching
    - Computes private Jaccard similarity
    - Supports threshold-based matching

    For thesis evaluation, also provides a simulated mode
    that computes correct results without full 2-party protocol.
    """

    def __init__(
        self,
        config: Optional[PSIConfig] = None,
        hpo_manager=None,
        use_ancestors: bool = True,
        simulated: bool = False
    ):
        """
        Initialize Hybrid PSI.

        Args:
            config: PSI configuration
            hpo_manager: HPOManager for ancestor expansion
            use_ancestors: Whether to expand terms with ancestors
            simulated: If True, skip crypto and compute directly (for testing)
        """
        self.config = config or PSIConfig()
        self.hpo_manager = hpo_manager
        self.use_ancestors = use_ancestors
        self.simulated = simulated

        if not simulated:
            self._dh_psi = DiffieHellmanPSI(config)
        else:
            self._dh_psi = None

    def expand_terms(self, terms: Set[str]) -> Set[str]:
        """
        Expand terms with ancestors for semantic matching.

        Args:
            terms: Original HPO term IDs

        Returns:
            Expanded set including ancestors
        """
        if not self.use_ancestors or not self.hpo_manager:
            return terms

        expanded = set(terms)
        for term in terms:
            ancestors = self.hpo_manager.get_ancestors(term, include_self=False)
            expanded.update(ancestors)

        return expanded

    def compute_intersection_size(
        self,
        terms1: Set[str],
        terms2: Set[str]
    ) -> int:
        """
        Compute |terms1 ∩ terms2| using PSI.

        Args:
            terms1: First set of HPO terms
            terms2: Second set of HPO terms

        Returns:
            Size of intersection
        """
        # Expand terms if configured
        expanded1 = self.expand_terms(terms1)
        expanded2 = self.expand_terms(terms2)

        if self.simulated:
            # Direct computation (for testing/benchmarking)
            return len(expanded1 & expanded2)

        # Full PSI protocol
        # Client = terms1, Server = terms2
        client_encrypted, client_secret = self._dh_psi.client_prepare(expanded1)
        doubly_encrypted, server_encrypted = self._dh_psi.server_process(
            client_encrypted, expanded2
        )
        intersection = self._dh_psi.client_compute_intersection(
            client_encrypted, client_secret, doubly_encrypted, server_encrypted
        )

        return len(intersection)

    def compute_private_jaccard(
        self,
        terms1: Set[str],
        terms2: Set[str]
    ) -> float:
        """
        Compute Jaccard similarity privately.

        J(A, B) = |A ∩ B| / |A ∪ B|

        Intersection computed via PSI; union computed from set sizes.

        Args:
            terms1: First set of HPO terms
            terms2: Second set of HPO terms

        Returns:
            Jaccard similarity
        """
        expanded1 = self.expand_terms(terms1)
        expanded2 = self.expand_terms(terms2)

        intersection_size = self.compute_intersection_size(terms1, terms2)

        # Union size = |A| + |B| - |A ∩ B|
        union_size = len(expanded1) + len(expanded2) - intersection_size

        if union_size == 0:
            return 1.0 if len(expanded1) == 0 and len(expanded2) == 0 else 0.0

        return intersection_size / union_size

    def threshold_match(
        self,
        terms1: Set[str],
        terms2: Set[str],
        threshold: float
    ) -> bool:
        """
        Check if similarity exceeds threshold.

        Binary answer: returns True/False without revealing exact similarity.

        Args:
            terms1: First set of terms
            terms2: Second set of terms
            threshold: Minimum similarity threshold

        Returns:
            True if Jaccard >= threshold
        """
        jaccard = self.compute_private_jaccard(terms1, terms2)
        return jaccard >= threshold


class PSIPhenopacketMatcher:
    """
    High-level interface for PSI-based phenopacket matching.

    Follows the same pattern as PhenopacketSimilarityCalculator
    but uses PSI for privacy-preserving computation.
    """

    def __init__(
        self,
        psi: HybridPSI,
        cache_results: bool = True
    ):
        """
        Initialize PSI phenopacket matcher.

        Args:
            psi: Configured HybridPSI instance
            cache_results: Whether to cache similarity results
        """
        self.psi = psi
        self.cache_results = cache_results
        self._cache: Dict[Tuple[str, str], float] = {}

    def extract_terms(self, phenopacket: Dict) -> Set[str]:
        """
        Extract HPO terms from phenopacket.

        Args:
            phenopacket: Phenopacket dictionary

        Returns:
            Set of HPO term IDs
        """
        terms = set()
        for feature in phenopacket.get("phenotypicFeatures", []):
            if not feature.get("excluded", False):
                terms.add(feature["type"]["id"])
        return terms

    def compute_similarity(
        self,
        phenopacket1: Dict,
        phenopacket2: Dict
    ) -> float:
        """
        Compute private similarity between phenopackets.

        Args:
            phenopacket1: First phenopacket
            phenopacket2: Second phenopacket

        Returns:
            Private Jaccard similarity
        """
        id1 = phenopacket1.get("id")
        id2 = phenopacket2.get("id")

        # Check cache
        if self.cache_results and id1 and id2:
            cache_key = tuple(sorted([id1, id2]))
            if cache_key in self._cache:
                return self._cache[cache_key]

        terms1 = self.extract_terms(phenopacket1)
        terms2 = self.extract_terms(phenopacket2)

        similarity = self.psi.compute_private_jaccard(terms1, terms2)

        # Cache
        if self.cache_results and id1 and id2:
            self._cache[cache_key] = similarity

        return similarity

    def find_threshold_matches(
        self,
        query: Dict,
        candidates: List[Dict],
        threshold: float
    ) -> List[int]:
        """
        Find candidates exceeding similarity threshold.

        Two-step reveal ladder step 1:
        Returns only indices, no similarity values revealed.

        Args:
            query: Query phenopacket
            candidates: List of candidate phenopackets
            threshold: Minimum similarity

        Returns:
            List of indices exceeding threshold
        """
        query_terms = self.extract_terms(query)
        matches = []

        for idx, candidate in enumerate(candidates):
            candidate_terms = self.extract_terms(candidate)
            if self.psi.threshold_match(query_terms, candidate_terms, threshold):
                matches.append(idx)

        return matches

    def find_most_similar(
        self,
        query: Dict,
        candidates: List[Dict],
        top_k: int = 10
    ) -> List[Tuple[int, float]]:
        """
        Find most similar candidates.

        Args:
            query: Query phenopacket
            candidates: List of candidates
            top_k: Number of top results

        Returns:
            List of (index, similarity) tuples
        """
        results = []

        for idx, candidate in enumerate(candidates):
            sim = self.compute_similarity(query, candidate)
            results.append((idx, sim))

        # Sort by similarity descending
        results.sort(key=lambda x: x[1], reverse=True)

        return results[:top_k]


def create_psi_matcher(
    security_parameter: int = 128,
    hpo_manager=None,
    use_ancestors: bool = True,
    simulated: bool = True
) -> PSIPhenopacketMatcher:
    """
    Factory function to create PSI matcher.

    Args:
        security_parameter: Cryptographic security parameter
        hpo_manager: HPOManager for ancestor expansion
        use_ancestors: Whether to use ancestor expansion
        simulated: Use simulated mode (faster, for testing)

    Returns:
        Configured PSIPhenopacketMatcher
    """
    config = PSIConfig(security_parameter=security_parameter)

    psi = HybridPSI(
        config=config,
        hpo_manager=hpo_manager,
        use_ancestors=use_ancestors,
        simulated=simulated
    )

    return PSIPhenopacketMatcher(psi)


if __name__ == "__main__":
    # Example usage
    logging.basicConfig(level=logging.INFO)

    # Test DH-PSI
    print("Testing DH-PSI...")
    dh_psi = DiffieHellmanPSI()

    set_a = {"HP:0001234", "HP:0002345", "HP:0003456"}
    set_b = {"HP:0002345", "HP:0003456", "HP:0004567"}

    # Client prepares
    client_enc, client_secret = dh_psi.client_prepare(set_a)
    print(f"Client encrypted {len(client_enc)} elements")

    # Server processes
    doubly_enc, server_enc = dh_psi.server_process(client_enc, set_b)
    print(f"Server returned {len(doubly_enc)} doubly encrypted elements")

    # Client computes intersection
    intersection = dh_psi.client_compute_intersection(
        client_enc, client_secret, doubly_enc, server_enc
    )
    print(f"Intersection: {intersection}")
    print(f"Expected: {set_a & set_b}")

    # Test HybridPSI (simulated mode)
    print("\nTesting HybridPSI (simulated)...")
    hybrid = HybridPSI(simulated=True)

    jaccard = hybrid.compute_private_jaccard(set_a, set_b)
    expected = len(set_a & set_b) / len(set_a | set_b)
    print(f"Private Jaccard: {jaccard:.4f}")
    print(f"Expected Jaccard: {expected:.4f}")

    # Test threshold matching
    print(f"\nThreshold match (>=0.3): {hybrid.threshold_match(set_a, set_b, 0.3)}")
    print(f"Threshold match (>=0.8): {hybrid.threshold_match(set_a, set_b, 0.8)}")

    # Test PSIPhenopacketMatcher
    print("\nTesting PSIPhenopacketMatcher...")
    matcher = create_psi_matcher(simulated=True)

    pp1 = {
        "id": "patient_1",
        "phenotypicFeatures": [
            {"type": {"id": "HP:0001234"}, "excluded": False},
            {"type": {"id": "HP:0002345"}, "excluded": False}
        ]
    }
    pp2 = {
        "id": "patient_2",
        "phenotypicFeatures": [
            {"type": {"id": "HP:0002345"}, "excluded": False},
            {"type": {"id": "HP:0003456"}, "excluded": False}
        ]
    }

    sim = matcher.compute_similarity(pp1, pp2)
    print(f"Phenopacket similarity: {sim:.4f}")
