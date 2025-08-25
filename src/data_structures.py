"""
Estruturas de dados tipadas para análise geométrica.
Substitui dicionários soltos por classes com validação e tipo.
"""

from dataclasses import dataclass, field
from typing import List, Dict, Optional, Tuple, Any
import math


@dataclass
class BoundingBox:
    """Representa um bounding box 3D."""
    xmin: float
    ymin: float
    zmin: float
    xmax: float
    ymax: float
    zmax: float
    
    @property
    def width(self) -> float:
        """Largura (eixo X)."""
        return abs(self.xmax - self.xmin)
    
    @property
    def height(self) -> float:
        """Altura (eixo Y)."""
        return abs(self.ymax - self.ymin)
    
    @property
    def depth(self) -> float:
        """Profundidade (eixo Z)."""
        return abs(self.zmax - self.zmin)
    
    @property
    def center(self) -> Tuple[float, float, float]:
        """Centro do bounding box."""
        return (
            (self.xmin + self.xmax) / 2,
            (self.ymin + self.ymax) / 2,
            (self.zmin + self.zmax) / 2
        )


@dataclass
class GeometricComponent:
    """CLASSE NÃO UTILIZADA - Mantida para compatibilidade. Representa um componente geométrico (cilindro, cone, etc.)."""
    geometric_type: str  # 'cilindrica', 'conica', 'plana'
    diameter: Optional[float] = None
    diameter_min: Optional[float] = None
    diameter_max: Optional[float] = None
    altura: float = 0.0
    area: float = 0.0
    direction: Optional[str] = None  # 'expanding', 'contracting'
    
    def __post_init__(self):
        """Validação após inicialização."""
        if self.geometric_type not in ['cilindrica', 'conica', 'plana']:
            raise ValueError(f"Tipo geométrico inválido: {self.geometric_type}")
        
        if self.geometric_type == 'conica':
            if self.diameter_min is None or self.diameter_max is None:
                raise ValueError("Cone deve ter diameter_min e diameter_max")
        elif self.geometric_type == 'cilindrica':
            if self.diameter is None:
                raise ValueError("Cilindro deve ter diameter")


@dataclass
class CircularFeature:
    """Representa um furo circular completo."""
    center: Tuple[float, float, float]
    min_d: float
    max_d: float
    has_cone: bool = False
    geometric_components: List[GeometricComponent] = field(default_factory=list)
    total_area: float = 0.0
    coordinates: Dict[str, Any] = field(default_factory=dict)
    
    @property
    def num_components(self) -> int:
        """Número de componentes geométricos."""
        return len(self.geometric_components)
    
    @property
    def is_simple(self) -> bool:
        """Verifica se é um furo simples (sem componentes complexos)."""
        return len(self.geometric_components) <= 1


@dataclass
class SemicircularFeature:
    """Representa um furo semicircular."""
    center: Tuple[float, float]
    min_d: float
    max_d: float
    group: Dict[str, Any] = field(default_factory=dict)
    geometric_components: Optional[List[GeometricComponent]] = None


@dataclass
class RectangularFeature:
    """Representa um furo retangular."""
    grupo_id: int
    num_faces: int
    comprimento: float
    largura: float
    faces: List[Any] = field(default_factory=list)
    bbox: Optional[BoundingBox] = None
    
    @property
    def area(self) -> float:
        """Área do furo retangular."""
        return self.comprimento * self.largura


@dataclass
class AnalysisWarning:
    """Representa um warning da análise."""
    warning_type: str  # 'circular', 'semicircular', 'rectangular'
    message: str
    count: int = 1


@dataclass
class PieceAnalysis:
    """CLASSE NÃO UTILIZADA - Mantida para compatibilidade. Resultado completo da análise de uma peça."""
    mold_name: str
    part_name: str
    dimensions: BoundingBox
    positions: int
    circular_features: List[CircularFeature] = field(default_factory=list)
    semicircular_features: List[SemicircularFeature] = field(default_factory=list)
    rectangular_features: List[RectangularFeature] = field(default_factory=list)
    warnings: List[AnalysisWarning] = field(default_factory=list)
    
    @property
    def total_circular_holes(self) -> int:
        """Total de furos circulares."""
        return len(self.circular_features)
    
    @property
    def total_semicircular_holes(self) -> int:
        """Total de furos semicirculares."""
        return len(self.semicircular_features)
    
    @property
    def total_rectangular_holes(self) -> int:
        """Total de furos retangulares."""
        return len(self.rectangular_features)


# Classe para tolerâncias e regras configuráveis
@dataclass
class AnalysisRules:
    """Regras e tolerâncias configuráveis para análise."""
    # Tolerâncias geométricas
    tol_face_border: float = 0.5
    tol_bbox: float = 1.0
    tol_planar_height: float = 1.0
    tol_conexao: float = 1e-3
    tol_center: float = 1.0
    tol_d: float = 0.5
    group_tol: float = 3.0
    lateral_tol: float = 150.0
    
    # Regras de profundidade mínima
    min_hole_depth: float = 1.0
    
    # Critérios para detecção de passagem
    passant_threshold_ratio: float = 0.95  # 95% da altura total
    
    # Tolerâncias de agrupamento
    center_grouping_tol: float = 15.0
    diameter_grouping_tol: float = 0.5


# Instância global das regras (pode ser configurada externamente)
DEFAULT_RULES = AnalysisRules()


def normalize_float(value: float, precision: int = 1) -> float:
    """
    FUNÇÃO NÃO UTILIZADA - Mantida para compatibilidade.
    Normaliza um valor float para agrupamento consistente.
    Encapsula a lógica de arredondamento para evitar inconsistências.
    """
    return round(value, precision)


def create_grouping_key(center: Tuple[float, float, float], diameter: float, 
                       rules: AnalysisRules = DEFAULT_RULES) -> Tuple[int, int, int]:
    """
    FUNÇÃO NÃO UTILIZADA - Mantida para compatibilidade.
    Cria chave de agrupamento consistente para features similares.
    """
    return (
        round(center[0] / rules.center_grouping_tol),
        round(center[1] / rules.center_grouping_tol),
        round(diameter / rules.diameter_grouping_tol)
    )
